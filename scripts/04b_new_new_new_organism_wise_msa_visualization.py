#!/usr/bin/env python3
"""
STAGE 04b — Organism-wise (genus-level) conservation profile + zoomed MSA
heatmap, built alongside the existing combined all-sequences visualization.

WHAT THIS DOES
--------------
For each genus 03b_new_organism_wise_phylogenetic_trees.py already aligned
(outputs/phylo_by_organism/<Genus>/<Genus>_aligned.fasta), this:
  1. Loads that genus's alignment directly -- NOT re-aligned here. 03b
     already ran a fresh per-genus MAFFT alignment for exactly this
     reason (a subset's optimal alignment differs from slicing rows out
     of the global MSA), so 04b consumes that file as-is.
  2. Runs the exact same conservation-profile + zoom-heatmap pipeline as
     stage 04, writing outputs into that SAME genus folder.
Then it ALSO (re-)builds the combined all-sequences visualization, using
stage 04's own ALIGNED_PATH/output paths, so running this ONE script
produces both the per-genus visualizations and the complete combined one,
mirroring 03b's design. Running stage 04 directly still works unchanged
and produces identical output -- nothing about stage 04 itself was
modified.

NEW, ADDITIVE: this ALSO now builds SPECIES-WISE (per-organism)
conservation profiles -- outputs/msa_by_organism/<Species>/. Unlike the
genus case, no stage in this pipeline produces a per-organism alignment
for this to reuse (03b only ever builds genus-pooled + one combined
alignment), so this writes and MAFFT-aligns each species' own sequence
subset itself, from outputs/family_sequences.fasta + outputs/
family_labels.csv's "organism" column -- same species-splitting
convention 03c's OrthoFinder step already uses. This is purely additive:
the genus-wise and combined outputs above are produced exactly as before,
byte-for-byte unchanged.

WHY NOTHING IS DUPLICATED HERE (unlike 03b)
--------------------------------------------
03b had to copy three functions verbatim from stage 03 because they get
pickled to multiprocessing worker processes under "spawn", and a module
loaded via importlib.util.spec_from_file_location() (required for
digit-prefixed filenames like "04_*.py") gets a synthetic module name
that a spawned child process can't re-import to unpickle against.

Stage 04 has no multiprocessing at all -- every function in it
(alignment_to_byte_array, compute_raw_conservation, smooth_conservation,
compute_gap_frequency, auto_select_zoom_window, plot_conservation_profile,
plot_zoom_heatmap) only ever runs synchronously in this same process. So
every one of them is imported directly from stage 04 below, with zero
duplication -- the same way 03b imports plot_and_save_nj_tree,
load_gene_labels, and clean_gene_name from stage 03 without copying them.

Requirements:
    pip install biopython numpy matplotlib
"""

import os
import re
import csv
import sys
import glob
import shutil
import subprocess
import importlib.util

from Bio import AlignIO, SeqIO

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE04_FILENAME = None   # set explicitly to skip auto-detection, e.g.
                           # "04_Multiple_sequence_alignment_visualization_conservation_profiling.py"

OUT_BASE_DIR = "outputs/phylo_by_organism"   # same base dir 03b writes into
MIN_SEQS_FOR_VISUALIZATION = 2               # a conservation profile needs at least 2 sequences;
                                              # in practice 03b already gates this more strictly
                                              # (MIN_SEQS_PER_GENUS) since an aligned file only
                                              # exists here if 03b produced one

FAMILY_FASTA_PATH = "outputs/family_sequences.fasta"   # stage 01's protein output -- MAFFT's input

# --------------------------------------------------------------------------
# NEW, additive: species-wise (per-organism) profiles. Unlike the genus
# case above, no stage in this pipeline produces a per-organism alignment
# anywhere (03b only ever builds genus-pooled + one combined alignment --
# see that script's own docstring), so this writes and aligns its own
# per-species FASTA subset, from scratch, into its own dedicated folder.
# Named "msa_by_organism" deliberately: this is also the exact folder-name
# substring the dashboard's "MSA Conservation Profile (species-wise)"
# panel pattern already searches for (*msa*by_organism*), so no dashboard
# change is needed once this runs.
# --------------------------------------------------------------------------
OUT_BASE_DIR_ORGANISM = "outputs/msa_by_organism"
LABELS_PATH = "outputs/family_labels.csv"        # confirmed "seq_id"/"organism" columns (stage 01)
MIN_SEQS_PER_ORGANISM = MIN_SEQS_FOR_VISUALIZATION   # same floor as genus/combined -- a
                                                       # conservation profile needs >= 2 sequences
                                                       # regardless of which grouping produced them
# ===========================================================================


def ensure_family_alignment(unaligned_path=FAMILY_FASTA_PATH, aligned_path=None):
    """If outputs/family_aligned.fasta is missing or empty, builds it with
    MAFFT directly -- rather than just skipping with instructions, per your
    request not to lose the combined-visualization figure.

    Self-contained on purpose (doesn't shell out to a separate
    02_align_with_mafft_for_all.py file): that script isn't guaranteed to
    exist next to this one, and this mirrors exactly what it does --
    `mafft --auto --quiet` -- so there's nothing to duplicate incorrectly.

    Includes the same WSL2 safeguard the original 02 script used: MAFFT
    writes working files to $TMPDIR (default /tmp), which on WSL2 is a
    small RAM-backed tmpfs that fills up fast and fails with a cryptic
    "No space left on device" even when the real disk has plenty of room.
    Pointed at a folder next to outputs/ instead, which lives on the real
    filesystem.

    Returns True if aligned_path exists and is non-empty afterward
    (whether it was already there or just got built), False otherwise --
    callers should treat False as "still can't proceed", not silently
    continue.
    """
    if aligned_path is None:
        aligned_path = msa_base.ALIGNED_PATH

    if os.path.exists(aligned_path) and os.path.getsize(aligned_path) > 0:
        return True

    print(f"[*] {aligned_path} not found or empty -- attempting to build it with MAFFT...")

    if shutil.which("mafft") is None:
        print("    [!] 'mafft' not found on PATH -- can't auto-build the alignment. "
              "Install it (conda install -c bioconda mafft) or run "
              "02_align_with_mafft_for_all.py yourself once it's installed.")
        return False

    if not (os.path.exists(unaligned_path) and os.path.getsize(unaligned_path) > 0):
        print(f"    [!] {unaligned_path} not found or empty either -- run the fetch "
              f"stage (01c/01d) first to produce it.")
        return False

    n_seqs = sum(1 for line in open(unaligned_path) if line.startswith(">"))
    print(f"    Aligning {n_seqs} sequences with MAFFT (this can take a while for "
          f"large families)...")

    tmp_dir = os.path.abspath("outputs/mafft_tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    env = os.environ.copy()
    env["TMPDIR"] = tmp_dir

    with open(aligned_path, "w") as out_f:
        result = subprocess.run(
            ["mafft", "--auto", "--quiet", unaligned_path],
            stdout=out_f, stderr=subprocess.PIPE, text=True, env=env,
        )

    if result.returncode != 0:
        print(f"    [!] MAFFT failed:\n{result.stderr}")
        if os.path.exists(aligned_path):
            os.remove(aligned_path)  # don't leave a stale/empty file behind
        return False

    if not (os.path.exists(aligned_path) and os.path.getsize(aligned_path) > 0):
        print(f"    [!] MAFFT reported success but {aligned_path} is missing or empty.")
        return False

    print(f"    [\u2713] Alignment built -> {aligned_path}")
    return True


def _find_stage04_module():
    if STAGE04_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE04_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE04_FILENAME is set to {STAGE04_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(glob.glob(os.path.join(_THIS_DIR, "04_*.py")))
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-04 visualization script (04_*.py) "
                  f"in {_THIS_DIR}. Set STAGE04_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 04_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE04_FILENAME above to override.\n")
    return candidates[0]


_STAGE04_PATH = _find_stage04_module()
_spec = importlib.util.spec_from_file_location("msa_base", _STAGE04_PATH)
msa_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(msa_base)  # only defines functions/constants -- its __main__
                                     # guard keeps its own combined-visualization build
                                     # from running here as a side effect of the import.


# ============================================================================
# New: genus-alignment discovery + shared per-alignment visualization pipeline
# ============================================================================
def discover_genus_alignments(base_dir):
    """Finds every outputs/phylo_by_organism/<Genus>/<Genus>_aligned.fasta
    03b already produced, so 04b visualizes exactly the same per-genus
    groupings used for the phylogenetic trees -- no separate genus-grouping
    logic duplicated here."""
    results = {}
    if not os.path.isdir(base_dir):
        return results
    for genus_dir in sorted(glob.glob(os.path.join(base_dir, "*"))):
        if not os.path.isdir(genus_dir):
            continue
        genus = os.path.basename(genus_dir)
        aligned_path = os.path.join(genus_dir, f"{genus}_aligned.fasta")
        if os.path.exists(aligned_path) and os.path.getsize(aligned_path) > 0:
            results[genus] = aligned_path
    return results


def visualize_alignment(aligned_path, out_profile_path, out_zoom_path, label=""):
    """Runs stage 04's exact pipeline (byte-array conversion, conservation
    scoring, smoothing, gap frequency, zoom-window selection, both plots)
    against ONE alignment file, writing to the given output paths. Every
    step here calls straight into msa_base -- nothing is reimplemented."""
    prefix = f"[{label}] " if label else ""

    if not (os.path.exists(aligned_path) and os.path.getsize(aligned_path) > 0):
        print(f"{prefix}[!] {aligned_path} not found or empty -- skipping.")
        return False

    alignment = AlignIO.read(aligned_path, "fasta")
    num_seqs = len(alignment)
    if num_seqs < MIN_SEQS_FOR_VISUALIZATION:
        print(f"{prefix}[!] Only {num_seqs} sequence(s) -- need >= "
              f"{MIN_SEQS_FOR_VISUALIZATION} for a meaningful conservation profile. Skipping.")
        return False

    seq_bytes = msa_base.alignment_to_byte_array(alignment)
    align_len = seq_bytes.shape[1]
    print(f"{prefix}Alignment loaded: {num_seqs} sequences x {align_len} columns")

    print(f"{prefix}Computing per-column conservation (vectorized)...")
    raw_scores = msa_base.compute_raw_conservation(seq_bytes, num_seqs)

    print(f"{prefix}Smoothing with sliding window (size={msa_base.WINDOW_SIZE})...")
    smoothed_scores = msa_base.smooth_conservation(raw_scores, msa_base.WINDOW_SIZE)

    print(f"{prefix}Computing gap frequency track...")
    gap_freq = msa_base.compute_gap_frequency(seq_bytes, num_seqs)

    x_positions = msa_base.np.arange(1, align_len + 1)

    if msa_base.ZOOM_START is not None and msa_base.ZOOM_END is not None:
        zoom_range = (msa_base.ZOOM_START - 1, msa_base.ZOOM_END)
        print(f"{prefix}Using manual zoom window: {msa_base.ZOOM_START}-{msa_base.ZOOM_END}")
    else:
        zoom_range = msa_base.auto_select_zoom_window(smoothed_scores, msa_base.ZOOM_WIDTH)
        print(f"{prefix}Auto-selected most conserved zoom window: "
              f"{zoom_range[0] + 1}-{zoom_range[1]}")

    msa_base.plot_conservation_profile(x_positions, smoothed_scores, gap_freq, num_seqs,
                                        align_len, zoom_range, out_profile_path)
    msa_base.plot_zoom_heatmap(alignment, seq_bytes, zoom_range, out_zoom_path)
    return True


def sanitize_organism_name(name):
    """Same alnum-only sanitizer used throughout this pipeline (01,
    03b/03c, 06b/06c, 08-organism, dashboard scripts) so folder names here
    match every other per-organism folder elsewhere."""
    if not name or str(name).strip() == "" or str(name).lower() == "nan":
        return "Unknown_organism"
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_")


def load_species_groups(labels_path):
    """Returns {sanitized_species_name: [seq_id, ...]} from
    family_labels.csv's 'organism' column -- same convention 03c uses for
    its own per-species OrthoFinder split."""
    if not os.path.exists(labels_path):
        return {}
    groups = {}
    with open(labels_path) as f:
        for row in csv.DictReader(f):
            seq_id = row.get("seq_id", "")
            organism = (row.get("organism") or "").strip()
            if not seq_id or not organism:
                continue
            species = sanitize_organism_name(organism)
            groups.setdefault(species, []).append(seq_id)
    return groups


def ensure_species_alignment(species, unaligned_path, aligned_path):
    """Per-species equivalent of ensure_family_alignment() above -- same
    MAFFT invocation, same flags, same WSL2 $TMPDIR safeguard -- kept as
    its own standalone copy rather than refactoring ensure_family_alignment
    to be reusable here, so that function's existing tested behavior stays
    completely untouched.

    Returns True if aligned_path exists and is non-empty afterward,
    False otherwise -- same contract as ensure_family_alignment()."""
    if os.path.exists(aligned_path) and os.path.getsize(aligned_path) > 0:
        return True

    if shutil.which("mafft") is None:
        print(f"    [!] {species}: 'mafft' not found on PATH -- can't build its alignment. "
              f"Install it: conda install -c bioconda mafft")
        return False

    if not (os.path.exists(unaligned_path) and os.path.getsize(unaligned_path) > 0):
        print(f"    [!] {species}: {unaligned_path} not found or empty -- nothing to align.")
        return False

    tmp_dir = os.path.abspath("outputs/mafft_tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    env = os.environ.copy()
    env["TMPDIR"] = tmp_dir

    with open(aligned_path, "w") as out_f:
        result = subprocess.run(
            ["mafft", "--auto", "--quiet", unaligned_path],
            stdout=out_f, stderr=subprocess.PIPE, text=True, env=env,
        )

    if result.returncode != 0:
        print(f"    [!] {species}: MAFFT failed:\n{result.stderr}")
        if os.path.exists(aligned_path):
            os.remove(aligned_path)
        return False

    if not (os.path.exists(aligned_path) and os.path.getsize(aligned_path) > 0):
        print(f"    [!] {species}: MAFFT reported success but {aligned_path} is missing or empty.")
        return False

    return True


def discover_or_build_species_alignments(base_dir, family_fasta_path, labels_path):
    """Species-wise equivalent of discover_genus_alignments() above, but
    since nothing upstream already produces a per-species alignment (see
    module docstring addendum), this builds them itself: writes each
    species' own sequence subset, then MAFFT-aligns it, reusing an
    already-built alignment on a re-run instead of redoing it (same
    idempotent-rerun convention as ensure_family_alignment()).

    Returns {sanitized_species_name: aligned_path} -- only species that
    ended up with a usable alignment."""
    species_groups = load_species_groups(labels_path)
    if not species_groups:
        print(f"[!] {labels_path} not found or empty -- can't determine species groupings, "
              f"skipping species-wise conservation profiles.")
        return {}

    if not (os.path.exists(family_fasta_path) and os.path.getsize(family_fasta_path) > 0):
        print(f"[!] {family_fasta_path} not found or empty -- skipping species-wise "
              f"conservation profiles.")
        return {}

    seq_lookup = SeqIO.to_dict(SeqIO.parse(family_fasta_path, "fasta"))

    results = {}
    for species, seq_ids in sorted(species_groups.items()):
        present_ids = [sid for sid in seq_ids if sid in seq_lookup]
        if len(present_ids) < MIN_SEQS_PER_ORGANISM:
            continue

        species_dir = os.path.join(base_dir, species)
        os.makedirs(species_dir, exist_ok=True)
        unaligned_path = os.path.join(species_dir, f"{species}_sequences.fasta")
        aligned_path = os.path.join(species_dir, f"{species}_aligned.fasta")

        if not (os.path.exists(aligned_path) and os.path.getsize(aligned_path) > 0):
            with open(unaligned_path, "w") as f:
                for sid in present_ids:
                    f.write(f">{sid}\n{str(seq_lookup[sid].seq)}\n")
            print(f"[*] {species}: aligning {len(present_ids)} sequence(s) with MAFFT...")
            if not ensure_species_alignment(species, unaligned_path, aligned_path):
                continue

        results[species] = aligned_path

    return results


def build_combined_visualization():
    """(Re-)builds the exact same combined all-sequences visualization
    stage 04 produces on its own, from the same ALIGNED_PATH to the same
    output paths -- so running THIS script alone gives you both outputs.
    Running stage 04 directly still works unchanged and produces
    identical output. Now auto-builds the alignment with MAFFT first if
    it's missing, rather than just skipping (see ensure_family_alignment)."""
    if not ensure_family_alignment():
        print(f"[!] Could not obtain {msa_base.ALIGNED_PATH} -- skipping the combined "
              f"all-sequences visualization.")
        return

    print(f"\n[*] Building the combined all-sequences conservation visualization "
          f"from {msa_base.ALIGNED_PATH} ...")
    visualize_alignment(msa_base.ALIGNED_PATH, msa_base.OUT_PROFILE_PATH,
                         msa_base.OUT_ZOOM_PATH, label="combined")


def main():
    genus_alignments = discover_genus_alignments(OUT_BASE_DIR)
    print(f"Found {len(genus_alignments)} genus/genera alignment(s) under '{OUT_BASE_DIR}': "
          f"{sorted(genus_alignments.keys())}\n")

    if not genus_alignments:
        print(f"[!] No <Genus>_aligned.fasta files found under {OUT_BASE_DIR} -- "
              f"run 03b_new_organism_wise_phylogenetic_trees.py first to produce them.")

    n_built, n_skipped = 0, 0
    for genus, aligned_path in sorted(genus_alignments.items()):
        genus_dir = os.path.dirname(aligned_path)
        out_profile = os.path.join(genus_dir, f"{genus}_conservation_profile.png")
        out_zoom = os.path.join(genus_dir, f"{genus}_msa_conservation_zoom.png")

        print(f"\n{'=' * 60}\n[\u27a4] {genus}\n{'=' * 60}")
        ok = visualize_alignment(aligned_path, out_profile, out_zoom, label=genus)
        n_built += int(ok)
        n_skipped += int(not ok)

    print(f"\n{'=' * 60}")
    print(f"Genus-wise conservation visualizations complete: {n_built} built, {n_skipped} skipped.")
    print(f"Outputs under: {OUT_BASE_DIR}/<Genus>/<Genus>_conservation_profile.png "
          f"and <Genus>_msa_conservation_zoom.png")
    print("=" * 60)

    # --------------------------------------------------------------------
    # NEW, additive: species-wise (per-organism). See module docstring
    # addendum and discover_or_build_species_alignments() for why this
    # builds its own alignments rather than discovering existing ones.
    # --------------------------------------------------------------------
    species_alignments = discover_or_build_species_alignments(
        OUT_BASE_DIR_ORGANISM, FAMILY_FASTA_PATH, LABELS_PATH
    )
    print(f"\nFound/built {len(species_alignments)} species alignment(s) under "
          f"'{OUT_BASE_DIR_ORGANISM}': {sorted(species_alignments.keys())}\n")

    n_built_sp, n_skipped_sp = 0, 0
    for species, aligned_path in sorted(species_alignments.items()):
        species_dir = os.path.dirname(aligned_path)
        out_profile = os.path.join(species_dir, f"{species}_conservation_profile.png")
        out_zoom = os.path.join(species_dir, f"{species}_msa_conservation_zoom.png")

        print(f"\n{'=' * 60}\n[\u27a4] {species}\n{'=' * 60}")
        ok = visualize_alignment(aligned_path, out_profile, out_zoom, label=species)
        n_built_sp += int(ok)
        n_skipped_sp += int(not ok)

    print(f"\n{'=' * 60}")
    print(f"Species-wise conservation visualizations complete: {n_built_sp} built, "
          f"{n_skipped_sp} skipped.")
    print(f"Outputs under: {OUT_BASE_DIR_ORGANISM}/<Species>/<Species>_conservation_profile.png "
          f"and <Species>_msa_conservation_zoom.png")
    print("=" * 60)

    build_combined_visualization()

    print("\n[\u2713] Done.")


if __name__ == "__main__":
    main()
