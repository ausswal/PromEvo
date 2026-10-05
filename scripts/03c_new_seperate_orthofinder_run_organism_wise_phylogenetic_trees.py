#!/usr/bin/env python3
"""
STAGE 03b — Organism-wise (genus-level) Neighbor-Joining phylogenetic trees,
built alongside the existing combined all-sequences tree.

WHAT THIS DOES
---------------
For each genus present in outputs/family_labels.csv's "organism" column
(e.g. all "Babesia ..." entries pooled together, per your Babesia example
and your choice of genus-level grouping over exact-species grouping), this:
  1. Pulls that genus's UNALIGNED sequences from outputs/family_sequences.fasta
     (the pre-alignment FASTA -- confirmed path/column names directly from
     01d_fetch_gene_family_multi_org.py's write_outputs(), not guessed).
  2. Re-aligns just that subset with MAFFT (a fresh alignment for the
     subset, not a slice of the global MSA -- a subset's optimal alignment
     genuinely differs from slicing rows out of an alignment computed
     across every genus together).
  3. Builds an NJ tree + bootstrap support using the EXACT SAME method as
     stage 03 (see "WHY FUNCTIONS ARE DUPLICATED BELOW", it matters).
Then it ALSO (re-)builds the combined all-sequences tree, using stage 03's
own ALIGNED_PATH/output paths, so running this ONE script produces both the
per-genus trees and the complete combined tree, as asked. Running stage 03
directly still works unchanged and produces the identical combined-tree
output -- nothing about stage 03 itself was modified.

NEW: after all of the above, this ALSO builds an independent REFERENCE
SPECIES TREE via OrthoFinder's STAG algorithm (Emms & Kelly, Genome Biology
2018) -- see build_orthofinder_reference_species_tree() below for the full
rationale. This is purely additive: every existing output (gene-level,
genus-level, combined tree) is produced exactly as before regardless of
whether OrthoFinder is installed or this step succeeds.

WHY THE MAFFT INVOCATION MIGHT NOT MATCH YOUR PIPELINE EXACTLY
-----------------------------------------------------------------
I have not seen 02_align_with_mafft_for_all.py's source, only its filename.
MAFFT_FLAGS below defaults to a standard "--auto" invocation. If your 02
script uses different flags (e.g. --localpair --maxiterate 1000 for higher
accuracy), change MAFFT_FLAGS to match exactly -- otherwise the per-genus
trees would be built with a subtly different alignment method than the
combined tree, which undermines comparing them side by side.

WHY THREE FUNCTIONS ARE DUPLICATED FROM STAGE 03 BELOW, VERBATIM, INSTEAD
OF BEING IMPORTED LIKE EVERYTHING ELSE
----------------------------------------------------------------------------
Stage 03 uses multiprocessing with the "spawn" start method for its
bootstrap replicates. Spawn re-imports a worker function's module BY NAME
in each freshly started child process. This script (like every other script
in this pipeline that reuses a sibling stage) has to load digit-prefixed
filenames like "03_Build_...py" via importlib.util.spec_from_file_location,
since Python's normal `import` statement can't handle a module name
starting with a digit. That gives the loaded module a synthetic name (e.g.
"tree_base") that is NOT a real importable module on disk under that name
-- so a spawned child process can't re-import it to unpickle a reference to
a function defined in it.

This was tested empirically before writing this file, not assumed: loading
a digit-prefixed script this way and handing one of its functions to
multiprocessing.Pool under spawn fails with
    PicklingError: Can't pickle <function ...>: import of module '...' failed
Running the same function when it's defined directly in the script that IS
being executed as __main__ (which is what stage 03 does, and what
multiprocessing's spawn machinery has a dedicated resolution path for)
works fine -- confirmed the same way.

So: build_nj_tree, _single_bootstrap_worker, and run_parallel_bootstrap_nj
are copied verbatim from stage 03 into this file (real duplication, not a
reimplementation -- kept byte-identical in logic) specifically because they
get pickled to worker processes. plot_and_save_nj_tree, load_gene_labels,
and clean_gene_name are only ever called synchronously in this process, so
those ARE imported from stage 03 directly below, with no duplication.

Requirements:
    pip install biopython numpy matplotlib
    mafft on PATH
"""

import os
import re
import csv
import sys
import json
import glob
import shutil
import subprocess
import multiprocessing as mp
import importlib.util

import numpy as np
import matplotlib
# Force the non-interactive Agg backend BEFORE pyplot is imported anywhere
# in this process -- including inside the dynamically-loaded stage03 module
# below (matplotlib.use() is a process-global setting, so this one call
# covers both). Without it, matplotlib defaults to TkAgg on systems where
# tkinter is available, and Tk's internal state gets torn down on the wrong
# thread when the spawned bootstrap worker processes exit, producing
# "RuntimeError: main thread is not in main loop" / leaked semaphores /
# "Aborted". Agg is the right choice here regardless -- this script only
# ever saves PNGs, it never shows an interactive window.
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE03_FILENAME = None   # set explicitly to skip auto-detection, e.g.
                           # "03_Build_neighbor_joining_phylogenetic_tree_with_thousand_bootstraps.py"

FAMILY_FASTA_PATH = "outputs/family_sequences.fasta"   # unaligned -- confirmed path from 01d
LABELS_PATH = "outputs/family_labels.csv"               # confirmed "organism" column from 01d

OUT_BASE_DIR = "outputs/phylo_by_organism"
MIN_SEQS_PER_GENUS = 4      # NJ + any meaningful bootstrap topology needs at least a few taxa

NUM_BOOTSTRAPS_PER_GENUS = None
# None -> reuse whatever NUM_BOOTSTRAPS stage 03 is currently set to, for
# exact methodological consistency with the combined tree. Set an int here
# to run genus trees with a different bootstrap count instead (e.g. lower,
# for speed across many small genus groups).

MAFFT_FLAGS = ["--auto"]   # see module docstring -- confirm this matches
                           # 02_align_with_mafft_for_all.py's actual flags

# --------------------------------------------------------------------------
# OrthoFinder reference species tree (NEW, additive -- see module docstring
# addendum below). This builds a SEPARATE species-level tree via
# OrthoFinder's STAG algorithm, which infers a species tree directly from
# multi-copy gene family data instead of averaging distances or picking one
# paralog as a stand-in. Everything above (gene-level, genus-level, and the
# combined tree) is completely unaffected by this -- it's an independent
# addition, not a replacement.
# --------------------------------------------------------------------------
DATA_MANIFEST_PATH = "outputs/.pipeline_data_manifest.json"   # written by the fetch stage

ORTHOFINDER_INPUT_DIR = "outputs/orthofinder_input"
REFERENCE_SPECIES_TREE_DIR = "outputs/reference_species_tree"
REFERENCE_SPECIES_TREE_NEWICK = os.path.join(REFERENCE_SPECIES_TREE_DIR, "orthofinder_species_tree.nwk")
REFERENCE_SPECIES_TREE_PNG = os.path.join(REFERENCE_SPECIES_TREE_DIR, "orthofinder_species_tree.png")

MIN_SPECIES_FOR_ORTHOFINDER = 3   # need at least a few taxa for a tree to mean anything
ORTHOFINDER_THREADS = max(1, (os.cpu_count() or 2) - 1)
# ===========================================================================


def _find_stage03_module():
    if STAGE03_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE03_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE03_FILENAME is set to {STAGE03_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(glob.glob(os.path.join(_THIS_DIR, "03_*.py")))
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-03 tree-building script (03_*.py) "
                  f"in {_THIS_DIR}. Set STAGE03_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 03_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE03_FILENAME above to override.\n")
    return candidates[0]


_STAGE03_PATH = _find_stage03_module()
_spec = importlib.util.spec_from_file_location("tree_base", _STAGE03_PATH)
tree_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tree_base)  # only defines functions/constants -- its __main__
                                      # guard keeps its own combined-tree build from
                                      # running here as a side effect of the import.


# ============================================================================
# Duplicated verbatim from stage 03 -- see module docstring for why.
# ============================================================================
def build_nj_tree(alignment):
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, method="nj")
    return constructor.build_tree(alignment)


def _single_bootstrap_worker(args):
    seq_records_data, align_len, num_seqs, seed = args
    np.random.seed(seed)
    resampled_indices = np.random.choice(align_len, size=align_len, replace=True)
    resampled_records = []
    for seq_id, seq_str in seq_records_data:
        boot_seq = "".join([seq_str[i] for i in resampled_indices])
        resampled_records.append(SeqRecord(Seq(boot_seq), id=seq_id, description=""))
    resampled_align = MultipleSeqAlignment(resampled_records)
    boot_tree = build_nj_tree(resampled_align)
    boot_clades = set()
    for clade in boot_tree.get_nonterminals():
        leaves = tuple(sorted([leaf.name for leaf in clade.get_terminals()]))
        if 1 < len(leaves) < num_seqs:
            boot_clades.add(leaves)
    return boot_clades


def run_parallel_bootstrap_nj(alignment, num_bootstraps=1000, label=""):
    """Identical to stage 03's function; `label` is a purely cosmetic
    addition to progress prints when building many genus trees in one run."""
    prefix = f"[{label}] " if label else ""
    print(f"{prefix}[*] Building primary Neighbor-Joining (NJ) tree...")
    main_tree = build_nj_tree(alignment)

    for clade in main_tree.get_nonterminals():
        clade.bootstrap_count = 0

    align_len = alignment.get_alignment_length()
    num_seqs = len(alignment)
    seq_records_data = [(rec.id, str(rec.seq)) for rec in alignment]

    num_cpus = max(1, mp.cpu_count() - 1)
    print(f"{prefix}[*] Dispatching {num_bootstraps} bootstrap replicates across {num_cpus} CPU threads...")

    worker_args = [
        (seq_records_data, align_len, num_seqs, np.random.randint(0, 1_000_000_000))
        for _ in range(num_bootstraps)
    ]

    completed = 0
    with mp.Pool(processes=num_cpus) as pool:
        for boot_clades in pool.imap_unordered(_single_bootstrap_worker, worker_args, chunksize=10):
            completed += 1
            if completed % 100 == 0 or completed == num_bootstraps:
                print(f"{prefix}  \u27a4 Progress: {completed}/{num_bootstraps} replicates done...")
            for clade in main_tree.get_nonterminals():
                main_leaves = tuple(sorted([leaf.name for leaf in clade.get_terminals()]))
                if main_leaves in boot_clades:
                    clade.bootstrap_count += 1

    for clade in main_tree.get_nonterminals():
        confidence = (clade.bootstrap_count / num_bootstraps) * 100
        clade.confidence = round(confidence, 1)

    return main_tree


# ============================================================================
# New: genus grouping + per-genus alignment/tree pipeline
# ============================================================================
def load_genus_groups(labels_path):
    """Returns {genus: [seq_id, ...]} from family_labels.csv's 'organism'
    column (confirmed exact column name -- see module docstring)."""
    if not os.path.exists(labels_path):
        sys.exit(f"{labels_path} not found -- run the fetch stage (01c/01d) first.")
    groups = {}
    with open(labels_path) as f:
        for row in csv.DictReader(f):
            seq_id = row.get("seq_id", "")
            organism = (row.get("organism") or "").strip()
            if not seq_id or not organism:
                continue
            genus = organism.split()[0]
            groups.setdefault(genus, []).append(seq_id)
    return groups


def run_mafft(input_fasta, output_fasta, flags=MAFFT_FLAGS):
    cmd = ["mafft"] + list(flags) + [input_fasta]
    print(f"    [*] Running MAFFT: {' '.join(cmd)} > {output_fasta}")
    with open(output_fasta, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0 or os.path.getsize(output_fasta) == 0:
        print(f"    [!] MAFFT failed for {input_fasta}. Stderr:\n{result.stderr}")
        return False
    return True


def clear_internal_names_for_confidence_output(tree):
    """Standard Newick can't cleanly carry both a node LABEL and a numeric
    support value at the same position -- Biopython's writer, when both are
    set, concatenates them into one unparseable string (verified: a clade
    with name="Inner1" and confidence=77.0 writes as "Inner177.00", which
    reads back as confidence=None). DistanceTreeConstructor auto-assigns
    "InnerN" names that aren't biologically meaningful, so -- following
    standard phylogenetics practice -- those are cleared here, right before
    writing, so ONLY the bootstrap support survives in the saved .nwk file.
    Mutates the tree in place; call this AFTER Phylo.draw()/PNG rendering
    (which reads confidence from the in-memory tree directly and is
    unaffected either way) and immediately before Phylo.write()."""
    for clade in tree.get_nonterminals():
        clade.name = None


def build_genus_tree(genus, seq_ids, seq_lookup, label_mapping, out_base_dir, num_bootstraps):
    missing = [sid for sid in seq_ids if sid not in seq_lookup]
    if missing:
        print(f"    [!] {len(missing)} seq_id(s) for {genus} listed in {LABELS_PATH} were not "
              f"found in {FAMILY_FASTA_PATH} (e.g. {missing[:3]}) -- proceeding with the rest.")
    present_ids = [sid for sid in seq_ids if sid in seq_lookup]

    if len(present_ids) < MIN_SEQS_PER_GENUS:
        print(f"[!] Skipping {genus}: only {len(present_ids)} sequence(s) available "
              f"(need >= {MIN_SEQS_PER_GENUS}).")
        return False

    genus_dir = os.path.join(out_base_dir, genus)
    os.makedirs(genus_dir, exist_ok=True)
    raw_fasta = os.path.join(genus_dir, f"{genus}_sequences.fasta")
    aligned_fasta = os.path.join(genus_dir, f"{genus}_aligned.fasta")
    out_newick = os.path.join(genus_dir, f"{genus}_nj_tree.nwk")
    out_png = os.path.join(genus_dir, f"{genus}_nj_tree.png")

    with open(raw_fasta, "w") as f:
        for sid in present_ids:
            f.write(f">{sid}\n{str(seq_lookup[sid].seq)}\n")

    print(f"\n{'=' * 60}\n[\u27a4] {genus}: {len(present_ids)} sequence(s)\n{'=' * 60}")
    if not run_mafft(raw_fasta, aligned_fasta):
        return False

    alignment = AlignIO.read(aligned_fasta, "fasta")
    nj_tree = run_parallel_bootstrap_nj(alignment, num_bootstraps=num_bootstraps, label=genus)

    clear_internal_names_for_confidence_output(nj_tree)
    Phylo.write(nj_tree, out_newick, "newick")
    print(f"    [\u2713] Newick tree -> {out_newick}")

    # plot_and_save_nj_tree (imported from stage 03) reads the MODULE-LEVEL
    # NUM_BOOTSTRAPS global for its title text, not a function argument --
    # set here so the title matches whatever count was ACTUALLY used for
    # this genus (safe: this function is only ever called synchronously in
    # this process, never pickled to a worker).
    tree_base.NUM_BOOTSTRAPS = num_bootstraps
    tree_base.plot_and_save_nj_tree(nj_tree, label_mapping, out_png)

    return True


def build_combined_tree(label_mapping):
    """(Re-)builds the exact same combined all-sequences tree stage 03
    produces on its own, from the same ALIGNED_PATH to the same output
    paths -- so running THIS script alone gives you both outputs. Running
    stage 03 directly still works unchanged and produces identical output."""
    if not (os.path.exists(tree_base.ALIGNED_PATH) and os.path.getsize(tree_base.ALIGNED_PATH) > 0):
        print(f"\n[!] {tree_base.ALIGNED_PATH} not found or empty -- skipping the combined "
              f"all-sequences tree. Run 02_align_with_mafft_for_all.py first.")
        return

    print(f"\n[*] Building the combined all-sequences tree from {tree_base.ALIGNED_PATH} ...")
    combined_alignment = AlignIO.read(tree_base.ALIGNED_PATH, "fasta")
    combined_num_bootstraps = tree_base.NUM_BOOTSTRAPS
    combined_tree = run_parallel_bootstrap_nj(
        combined_alignment, num_bootstraps=combined_num_bootstraps, label="combined")

    clear_internal_names_for_confidence_output(combined_tree)
    Phylo.write(combined_tree, tree_base.OUT_NEWICK_TREE, "newick")
    print(f"[\u2713] Combined Newick tree -> {tree_base.OUT_NEWICK_TREE}")

    tree_base.NUM_BOOTSTRAPS = combined_num_bootstraps
    tree_base.plot_and_save_nj_tree(combined_tree, label_mapping, tree_base.OUT_TREE_IMAGE)


# ============================================================================
# NEW: OrthoFinder reference species tree
#
# WHY THIS IS ITS OWN INDEPENDENT STEP, NOT A MODIFICATION OF ANYTHING ABOVE
# ----------------------------------------------------------------------------
# Everything above this point builds trees where every SEQUENCE (every
# paralog) is its own tip -- that's correct for looking at gene-level
# relationships, but doesn't answer "how do the host species relate to each
# other", since a species with 5 paralogs shows up as 5 separate tips
# scattered around the tree rather than one coherent position.
#
# OrthoFinder's STAG algorithm (Emms & Kelly, Genome Biology 2018) exists
# specifically to build a SPECIES tree from data where genes have been
# duplicated -- it infers orthogroups (clusters of orthologous + paralogous
# sequences) and builds the species tree from the signal across all of them,
# rather than requiring one single-copy representative per species. That's
# a meaningfully different (and more principled) answer than either picking
# one paralog as a stand-in, or naively averaging pairwise distances.
#
# This pipeline only fetches ONE gene family (not whole proteomes), so this
# runs OrthoFinder on a per-species split of THAT family's protein
# sequences -- not a genome-wide analysis. With few species or very few
# paralogs per species this will be a weakly-supported tree (that's a
# real data limitation, not a bug); OrthoFinder will still run and produce
# its best estimate.
#
# REQUIRES REAL PROTEIN SEQUENCES. OrthoFinder clusters by protein
# similarity (Diamond blastp internally) -- running it on nucleotide gene
# sequences would silently produce meaningless orthogroups (nucleotide
# k-mers aren't what its similarity search expects). This checks
# DATA_MANIFEST_PATH (written by the fetch stage's protein/nucleotide
# detection -- see 01e_new5_fetch_family_and_promoters.py's
# write_data_manifest()) before running anything.
# ============================================================================

def sanitize_organism_name(name):
    """Same alnum-only sanitizer used throughout the pipeline (01, 06b/06c,
    08-organism, 011h dashboard) so this matches folder/label conventions
    used everywhere else."""
    if not name or str(name).strip() == "" or str(name).lower() == "nan":
        return "Unknown_organism"
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_")


def load_species_groups_by_genus(labels_path):
    """Returns {genus: {sanitized_species_name: [seq_id, ...]}} -- the
    same data load_species_groups() produces, but partitioned by genus
    first, so OrthoFinder can be run SEPARATELY per genus.

    Why this matters (and why the combined run alone isn't enough): a
    single OrthoFinder run over Plasmodium + Toxoplasma + Babesia together
    infers orthogroups across all three at once. Those genera are distant
    enough that many orthogroups will contain sequences from only one
    genus, and the resulting STAG species tree's within-genus branching is
    inferred in the presence of (and partly constrained by) the deep
    between-genus splits. Running each genus separately means each tree's
    orthogroups, and therefore its within-genus topology, are inferred
    only from that genus's own proteins -- which is what you want when the
    question is "how do species relate WITHIN this genus", rather than
    "how do these genera relate to each other".

    Both are produced (see main()): per-genus trees for within-genus
    relationships, plus the original combined tree for the across-genus
    view. Neither replaces the other.
    """
    if not os.path.exists(labels_path):
        return {}
    by_genus = {}
    with open(labels_path) as f:
        for row in csv.DictReader(f):
            seq_id = row.get("seq_id", "")
            organism = (row.get("organism") or "").strip()
            if not seq_id or not organism:
                continue
            genus = organism.split()[0]
            species = sanitize_organism_name(organism)
            by_genus.setdefault(genus, {}).setdefault(species, []).append(seq_id)
    return by_genus


def load_species_groups(labels_path):
    """Returns {sanitized_species_name: [seq_id, ...]} from
    family_labels.csv's 'organism' column -- the FULL organism name this
    time (species-level), not just the genus, unlike load_genus_groups()."""
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


def check_real_protein_available():
    """Reads the fetch stage's data manifest. Defaults to True (proceed) if
    the manifest is missing -- same safe default used everywhere else this
    manifest is read in this pipeline (master_script_*.py's
    load_data_manifest()) -- so older runs from before the manifest existed
    aren't blocked."""
    if not os.path.exists(DATA_MANIFEST_PATH):
        print(f"[!] {DATA_MANIFEST_PATH} not found -- can't confirm these are real "
              f"protein sequences. Proceeding anyway (this pipeline's default "
              f"assumption when the manifest is missing); OrthoFinder results will "
              f"be meaningless if this data turns out to be nucleotide.")
        return True
    try:
        with open(DATA_MANIFEST_PATH) as f:
            manifest = json.load(f)
        return manifest.get("has_real_protein", True)
    except Exception as e:
        print(f"[!] Could not read {DATA_MANIFEST_PATH} ({e}) -- proceeding anyway.")
        return True


def write_orthofinder_species_fastas(species_groups, seq_lookup, input_dir):
    """Writes one protein FASTA per species (OrthoFinder's taxon name comes
    directly from the filename stem, so sanitized species names are used
    as filenames -- consistent with every other per-species folder/file
    name in this pipeline). Clears input_dir first so a species removed
    since the last run doesn't linger and get included by mistake."""
    if os.path.exists(input_dir):
        shutil.rmtree(input_dir)
    os.makedirs(input_dir, exist_ok=True)

    written = {}
    for species, seq_ids in sorted(species_groups.items()):
        present_ids = [sid for sid in seq_ids if sid in seq_lookup]
        if not present_ids:
            continue
        fasta_path = os.path.join(input_dir, f"{species}.fasta")
        with open(fasta_path, "w") as f:
            for sid in present_ids:
                f.write(f">{sid}\n{str(seq_lookup[sid].seq)}\n")
        written[species] = len(present_ids)

    return written


def find_latest_orthofinder_results(input_dir):
    """OrthoFinder writes into <input_dir>/OrthoFinder/Results_<date>/, with
    a folder name that changes every run (and gets a _1, _2... suffix if
    run more than once on the same day) -- this finds the most recently
    modified one rather than assuming a fixed name."""
    pattern = os.path.join(input_dir, "OrthoFinder", "Results_*")
    candidates = glob.glob(pattern)
    if not candidates:
        return None
    return max(candidates, key=os.path.getmtime)


def locate_species_tree_newick(results_dir):
    """OrthoFinder names this file consistently across versions; the
    node-labeled variant (support values named on internal nodes) is
    preferred when present."""
    candidates = [
        os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt"),
        os.path.join(results_dir, "Species_Tree", "SpeciesTree_rooted.txt"),
    ]
    for path in candidates:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            return path
    return None


def plot_reference_species_tree(newick_path, out_png, title_suffix=""):
    """Simple, self-contained rendering -- deliberately NOT reusing
    tree_base.plot_and_save_nj_tree(), since that function's title/caption
    is written specifically for stage 03's bootstrap-replicate NJ trees
    (it says "N bootstrap replicates"), which would be a factually wrong
    caption on an OrthoFinder/STAG tree that uses a completely different
    support-value method."""
    tree = Phylo.read(newick_path, "newick")
    for tip in tree.get_terminals():
        tip.name = tip.name.replace("_", " ") if tip.name else tip.name

    fig = plt.figure(figsize=(10, max(4, 0.35 * tree.count_terminals() + 2)))
    ax = fig.add_subplot(1, 1, 1)
    title = "Species Tree (OrthoFinder / STAG)"
    if title_suffix:
        title = f"{title_suffix} — {title}"
    ax.set_title(title, fontsize=13, fontweight="bold")
    Phylo.draw(tree, axes=ax, do_show=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def _run_orthofinder_for_group(species_groups, seq_lookup, input_dir, out_newick,
                                out_png, group_label, min_species=MIN_SPECIES_FOR_ORTHOFINDER):
    """Runs one complete OrthoFinder job for a given set of species and
    copies out its species tree. Shared by BOTH the combined run and each
    per-genus run (extracted rather than duplicated so the two can never
    drift apart in behavior). Returns True/False; never raises."""
    if len(species_groups) < 2:
        print(f"[!] {group_label}: only {len(species_groups)} species -- need at least 2 "
              f"to build any tree. Skipping.")
        return False
    if len(species_groups) < min_species:
        print(f"[!] {group_label}: only {len(species_groups)} species (recommended >= "
              f"{min_species}) -- proceeding anyway, but treat this tree's support "
              f"values with caution.")

    written = write_orthofinder_species_fastas(species_groups, seq_lookup, input_dir)
    if len(written) < 2:
        print(f"[!] {group_label}: only {len(written)} species had sequences actually "
              f"present in {FAMILY_FASTA_PATH} -- skipping.")
        return False
    print(f"[*] {group_label}: wrote {len(written)} per-species protein FASTA(s) to "
          f"{input_dir}: {written}")

    print(f"[*] {group_label}: running OrthoFinder (this can take a few minutes)...")
    result = subprocess.run(
        ["orthofinder", "-f", input_dir, "-t", str(ORTHOFINDER_THREADS),
         "-a", str(ORTHOFINDER_THREADS)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"[!] {group_label}: OrthoFinder exited with code {result.returncode}. "
              f"Stderr:\n{result.stderr[-2000:]}")
        return False

    results_dir = find_latest_orthofinder_results(input_dir)
    if not results_dir:
        print(f"[!] {group_label}: OrthoFinder ran but no Results_* folder was found "
              f"under {input_dir}/OrthoFinder/. Full stdout:\n{result.stdout[-2000:]}")
        return False

    newick_path = locate_species_tree_newick(results_dir)
    if not newick_path:
        print(f"[!] {group_label}: OrthoFinder completed but no species tree file was "
              f"found under {results_dir}/Species_Tree/. This can happen with very few "
              f"orthogroups -- check {results_dir} directly.")
        return False

    os.makedirs(os.path.dirname(out_newick), exist_ok=True)
    shutil.copy(newick_path, out_newick)
    print(f"[\u2713] {group_label}: species tree (Newick) -> {out_newick}")

    try:
        plot_reference_species_tree(newick_path=out_newick, out_png=out_png,
                                     title_suffix=group_label)
        print(f"[\u2713] {group_label}: species tree (image) -> {out_png}")
    except Exception as e:
        print(f"[!] {group_label}: tree built successfully but rendering the PNG failed: {e}")
        print(f"    The Newick file at {out_newick} is still valid and usable.")

    print(f"[*] {group_label}: full OrthoFinder results kept at: {results_dir}")
    return True


def build_orthofinder_per_genus_species_trees(seq_lookup):
    """NEW: runs OrthoFinder SEPARATELY for each genus, so each genus's
    species tree is inferred only from that genus's own proteins -- see
    load_species_groups_by_genus()'s docstring for why this differs
    meaningfully from the combined run, which is still performed too.

    Outputs go to outputs/reference_species_tree/by_genus/<Genus>/, keeping
    them clearly separate from the combined tree's own output location."""
    print(f"\n{'#' * 60}")
    print("# PER-GENUS SPECIES TREES (OrthoFinder / STAG)")
    print(f"{'#' * 60}")

    if shutil.which("orthofinder") is None:
        print("[!] 'orthofinder' not found on PATH -- skipping per-genus species trees. "
              "Install it with: conda install -c bioconda orthofinder")
        return

    if not check_real_protein_available():
        print("[!] Data is not confirmed to be real protein sequences -- skipping "
              "per-genus species trees (OrthoFinder requires protein input).")
        return

    by_genus = load_species_groups_by_genus(LABELS_PATH)
    if not by_genus:
        print(f"[!] No genus/species groups could be read from {LABELS_PATH}. Skipping.")
        return

    print(f"[*] {len(by_genus)} genus/genera found: "
          f"{ {g: len(sp) for g, sp in by_genus.items()} }")

    n_built, n_skipped = 0, 0
    for genus, species_groups in sorted(by_genus.items()):
        print(f"\n{'=' * 60}\n[\u27a4] Genus: {genus} ({len(species_groups)} species)\n{'=' * 60}")

        if len(species_groups) < 2:
            print(f"[i] {genus}: only {len(species_groups)} species in this genus -- a "
                  f"species tree needs at least 2 taxa, so there's nothing to infer "
                  f"here. Skipping (this is expected, not an error, when a genus is "
                  f"represented by a single species).")
            n_skipped += 1
            continue

        genus_out_dir = os.path.join(REFERENCE_SPECIES_TREE_DIR, "by_genus", genus)
        ok = _run_orthofinder_for_group(
            species_groups=species_groups,
            seq_lookup=seq_lookup,
            input_dir=os.path.join(ORTHOFINDER_INPUT_DIR + "_by_genus", genus),
            out_newick=os.path.join(genus_out_dir, f"{genus}_orthofinder_species_tree.nwk"),
            out_png=os.path.join(genus_out_dir, f"{genus}_orthofinder_species_tree.png"),
            group_label=genus,
        )
        n_built += int(ok)
        n_skipped += int(not ok)

    print(f"\n{'=' * 60}")
    print(f"Per-genus species trees complete: {n_built} built, {n_skipped} skipped.")
    print(f"Outputs under: {REFERENCE_SPECIES_TREE_DIR}/by_genus/<Genus>/")
    print("=" * 60)


def build_orthofinder_reference_species_tree(seq_lookup):
    """Top-level orchestration. Returns True/False; NEVER raises -- a
    failure here (OrthoFinder not installed, too few species, nucleotide
    data, etc.) is reported clearly and this simply gets skipped, leaving
    every other tree this script builds completely unaffected."""

    print(f"\n{'#' * 60}")
    print("# REFERENCE SPECIES TREE (OrthoFinder / STAG)")
    print(f"{'#' * 60}")

    if shutil.which("orthofinder") is None:
        print(
            "[!] 'orthofinder' not found on PATH -- skipping the reference species "
            "tree. Install it with: conda install -c bioconda orthofinder"
        )
        return False

    if not check_real_protein_available():
        print("[!] Data is not confirmed to be real protein sequences -- skipping "
              "the reference species tree (OrthoFinder requires protein input).")
        return False

    species_groups = load_species_groups(LABELS_PATH)
    return _run_orthofinder_for_group(
        species_groups=species_groups,
        seq_lookup=seq_lookup,
        input_dir=ORTHOFINDER_INPUT_DIR,
        out_newick=REFERENCE_SPECIES_TREE_NEWICK,
        out_png=REFERENCE_SPECIES_TREE_PNG,
        group_label="All genera combined",
    )


def main():
    mp.set_start_method("spawn", force=True)

    if not os.path.exists(FAMILY_FASTA_PATH):
        sys.exit(f"{FAMILY_FASTA_PATH} not found -- run the fetch stage (01c/01d) first.")

    label_mapping = tree_base.load_gene_labels(LABELS_PATH)
    seq_lookup = SeqIO.to_dict(SeqIO.parse(FAMILY_FASTA_PATH, "fasta"))
    genus_groups = load_genus_groups(LABELS_PATH)

    num_bootstraps = (NUM_BOOTSTRAPS_PER_GENUS if NUM_BOOTSTRAPS_PER_GENUS is not None
                       else tree_base.NUM_BOOTSTRAPS)

    counts = {g: len(ids) for g, ids in genus_groups.items()}
    print(f"Found {len(genus_groups)} genus/genera in {LABELS_PATH}: {counts}")
    print(f"Building genus-wise trees with {num_bootstraps} bootstrap replicate(s) each "
          f"(MIN_SEQS_PER_GENUS={MIN_SEQS_PER_GENUS}).\n")

    os.makedirs(OUT_BASE_DIR, exist_ok=True)
    n_built, n_skipped = 0, 0
    for genus, seq_ids in sorted(genus_groups.items()):
        ok = build_genus_tree(genus, seq_ids, seq_lookup, label_mapping, OUT_BASE_DIR, num_bootstraps)
        n_built += int(ok)
        n_skipped += int(not ok)

    print(f"\n{'=' * 60}")
    print(f"Genus-wise trees complete: {n_built} built, {n_skipped} skipped "
          f"(too few sequences or MAFFT/parse failure).")
    print(f"Outputs under: {OUT_BASE_DIR}/<Genus>/")
    print("=" * 60)

    build_combined_tree(label_mapping)

    # NEW, additive -- see function docstring. Never affects anything above:
    # failure here is reported and skipped, not raised.
    build_orthofinder_per_genus_species_trees(seq_lookup)
    build_orthofinder_reference_species_tree(seq_lookup)


if __name__ == "__main__":
    main()
