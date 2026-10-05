#!/usr/bin/env python3
"""
STAGE 08-organism — Organism-wise AND Genus-wise Promoter CpG Island
Discovery & Comparative Analysis.

WHY PER-GROUP COMPUTATION MATTERS HERE (a DIFFERENT reason than 06b/06c)
---------------------------------------------------------------------------
06b/06c's per-genus necessity is about a STOCHASTIC control (a shuffled
background) needing enough sequences to be statistically meaningful. That's
NOT the issue here -- 08's adaptive threshold selection (standard
Gardiner-Garden >=50% GC vs. the AT-rich fallback, and the AT-rich cutoff's
exact value) is a DETERMINISTIC function of the group's own global GC%, so
it's always well-defined regardless of group size.

The real problem with pooling everything (what 08 does today) is different:
computing ONE global GC% across every organism and picking ONE threshold
for the whole dataset means a 64%-GC organism and a 21%-GC organism (a real
spread observed in this pipeline, per the companion 06a/06b-organism
scripts) get forced through the identical cutoff. That defeats the entire
point of the adaptive-threshold design, which exists specifically because
different organisms sit at very different baseline GC content. Per-group
computation is what actually lets the adaptive logic do its job.

Because of that difference, MIN_SEQS_PER_GROUP defaults to 1 here (i.e. no
group is skipped by default) -- unlike 06b/06c/07c-organism, there's no
"too few sequences to be statistically valid" argument to justify skipping
a small group. Raise it yourself if you want to exclude groups too small to
be a meaningful UNIT OF COMPARISON (e.g. a single stray sequence dominating
one genus's whole plot), but that's a presentation choice, not a
statistical requirement, so it isn't imposed by default.

BOTH RESOLUTION LEVELS, ONE RUN
---------------------------------------------------------------------------
Unlike 06b/06c/07c-organism (genus-only), this produces BOTH:
  outputs/CpG_Island_By_Organism/<Organism>/   (full organism+strain, e.g.
                                                 Babesia_bovis_T2Bo)
  outputs/CpG_Island_By_Genus/<Genus>/         (e.g. Babesia)
Expect the two levels' min_gc_threshold to sometimes disagree for the same
species -- that's not a bug, it's the adaptive design responding correctly
to composition differences within a genus. A large organism-vs-genus
threshold gap for one species is itself a signal of within-genus GC
heterogeneity worth a second look.

REUSE
------
calculate_metrics, compute_global_nucleotide_content,
scan_promoter_for_cpg_islands, export_bed_file, and plot_cpg_summary are
imported directly from 08, unchanged. scan_promoter_for_cpg_islands takes a
FASTA path (not an in-memory sequence dict), so each group's sequences are
written to their own small FASTA file first -- same approach used in
07c-organism for MEME, for the same reason (the base function's signature
requires a real file).

plot_cpg_summary() writes to a MODULE-LEVEL GLOBAL (OUT_PLOT), not a
parameter -- exactly the same situation as 07c's TOMTOM_OUT_DIR. This
script patches cpb.OUT_PLOT to each group's own folder immediately before
calling cpb.plot_cpg_summary() for that group.

Requirements:
    pip install pandas numpy matplotlib seaborn biopython
"""

import os
import sys
import glob
import importlib.util

import pandas as pd
from Bio import SeqIO

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE08_FILENAME = None   # set explicitly to skip auto-detection

PROMOTER_FASTA_PATH = "outputs/promoter_sequences.fasta"   # same source 08 already uses
LABELS_PATH = "outputs/family_labels.csv"                   # confirmed "organism" column (01)

OUT_ROOT_ORGANISM = "outputs/CpG_Island_By_Organism"
OUT_ROOT_GENUS = "outputs/CpG_Island_By_Genus"

MIN_SEQS_PER_GROUP = 1   # see module docstring for why this differs from 06b/06c/07c-organism
# ===========================================================================


def _find_stage08_module():
    if STAGE08_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE08_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE08_FILENAME is set to {STAGE08_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(
        p for p in glob.glob(os.path.join(_THIS_DIR, "08_*.py"))
        if os.path.basename(p) != os.path.basename(__file__)
    )
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-08 CpG island script (08_*.py) "
                  f"in {_THIS_DIR}. Set STAGE08_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 08_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE08_FILENAME above to override.\n")
    return candidates[0]


_STAGE08_PATH = _find_stage08_module()
_spec = importlib.util.spec_from_file_location("cpg_base", _STAGE08_PATH)
cpb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cpb)  # only defines functions/constants -- its __main__
                                 # guard keeps its own pooled run from happening
                                 # here as an import side effect.

os.makedirs(OUT_ROOT_ORGANISM, exist_ok=True)
os.makedirs(OUT_ROOT_GENUS, exist_ok=True)


# ==============================================================================
# Dual-resolution data loader -- reads the FASTA + labels ONCE, groups by
# both organism and genus in the same pass.
# ==============================================================================
def sanitize_organism_name(name):
    """Same alnum-only sanitizer used throughout the pipeline (01, 05,
    23_comparative_matrix.py, 06b/06c/07c-organism) so labels/folder names
    line up everywhere."""
    import re
    if not name or str(name).strip() == "" or str(name).lower() == "nan":
        return "Unknown_organism"
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_")


def load_promoters_dual_grouped():
    """Returns (organism_groups, genus_groups), each {label: {seq_id: seq}}."""
    if not os.path.exists(PROMOTER_FASTA_PATH):
        sys.exit(f"[!] {PROMOTER_FASTA_PATH} not found -- run the promoter-fetch stage first.")
    if not os.path.exists(LABELS_PATH):
        sys.exit(f"[!] {LABELS_PATH} not found -- run the fetch stage (01) first "
                 f"(grouping needs its 'organism' column).")

    labels_df = pd.read_csv(LABELS_PATH)
    if "organism" not in labels_df.columns:
        sys.exit(f"[!] {LABELS_PATH} has no 'organism' column -- can't group by organism/genus. "
                 f"Actual columns: {list(labels_df.columns)}")
    organism_by_seqid = dict(zip(labels_df["seq_id"], labels_df["organism"]))

    all_records = list(SeqIO.parse(PROMOTER_FASTA_PATH, "fasta"))
    if not all_records:
        sys.exit(f"[!] No sequences parsed from {PROMOTER_FASTA_PATH}.")

    organism_groups, genus_groups = {}, {}
    n_unmapped = 0
    for rec in all_records:
        organism = organism_by_seqid.get(rec.id)
        if not organism or not str(organism).strip():
            n_unmapped += 1
            continue
        org_label = sanitize_organism_name(organism)
        genus_label = sanitize_organism_name(str(organism).strip().split()[0])
        seq = str(rec.seq)
        organism_groups.setdefault(org_label, {})[rec.id] = seq
        genus_groups.setdefault(genus_label, {})[rec.id] = seq

    if n_unmapped:
        print(f"[!] {n_unmapped} promoter sequence(s) in {PROMOTER_FASTA_PATH} had no "
              f"matching/non-empty 'organism' in {LABELS_PATH} and were excluded from both groupings.")

    for groups, label_word in ((organism_groups, "organism"), (genus_groups, "genus")):
        skipped = {g: len(s) for g, s in groups.items() if len(s) < MIN_SEQS_PER_GROUP}
        for g, n in skipped.items():
            print(f"[!] Skipping {label_word} '{g}': only {n} sequence(s) "
                  f"(< MIN_SEQS_PER_GROUP={MIN_SEQS_PER_GROUP}).")
            del groups[g]

    if not organism_groups and not genus_groups:
        sys.exit(f"[!] No organism or genus had >= {MIN_SEQS_PER_GROUP} promoter sequence(s).")
    return organism_groups, genus_groups


# ==============================================================================
# Per-group run -- a direct port of 08's own main() body, parameterized by
# group instead of hardcoded to the pooled dataset.
# ==============================================================================
def run_one_group(group_label, seqs, out_root, group_type_word):
    n_seqs = len(seqs)
    print(f"\n{'=' * 60}\n[{group_type_word}] {group_label}: {n_seqs} promoter sequence(s)\n{'=' * 60}")

    group_dir = os.path.join(out_root, group_label)
    os.makedirs(group_dir, exist_ok=True)

    group_fasta = os.path.join(group_dir, f"{group_label}_promoter_sequences.fasta")
    with open(group_fasta, "w") as f:
        for seq_id, seq in seqs.items():
            f.write(f">{seq_id}\n{seq}\n")

    df_islands, df_summary, min_gc_threshold = cpb.scan_promoter_for_cpg_islands(group_fasta)

    detailed_csv = os.path.join(group_dir, f"{group_label}_cpg_islands_detailed.csv")
    summary_csv = os.path.join(group_dir, f"{group_label}_cpg_island_summary_per_gene.csv")
    bed_path = os.path.join(group_dir, f"{group_label}_cpg_islands.bed")
    plot_path = os.path.join(group_dir, f"{group_label}_cpg_island_distribution.png")

    df_islands.to_csv(detailed_csv, index=False)
    df_summary.to_csv(summary_csv, index=False)
    cpb.export_bed_file(df_islands, bed_path)

    total_genes = len(df_summary)
    genes_with_islands = df_summary['has_cpg_island'].sum() if total_genes > 0 else 0
    print(f"    Total Promoter Sequences Analyzed : {total_genes}")
    print(f"    Promoters Containing CpG Island(s): {genes_with_islands} "
          f"({(genes_with_islands / total_genes * 100 if total_genes else 0):.1f}%)")
    print(f"    Total CpG Islands Identified      : {len(df_islands)}")
    print(f"    [\u2713] Detailed -> {detailed_csv}")
    print(f"    [\u2713] Summary  -> {summary_csv}")

    # Redirect the base module's OUT_PLOT global so plot_cpg_summary()
    # writes into THIS group's own folder (see module docstring) instead of
    # the pooled outputs/ location it defaults to.
    cpb.OUT_PLOT = plot_path
    cpb.plot_cpg_summary(df_summary, df_islands, min_gc_threshold)


def main():
    print("=" * 65)
    print("  ORGANISM-WISE + GENUS-WISE PROMOTER CpG ISLAND ANALYSIS  ")
    print("=" * 65)

    organism_groups, genus_groups = load_promoters_dual_grouped()

    print(f"\n{'#' * 60}\n# ORGANISM-LEVEL CpG ISLAND ANALYSIS\n{'#' * 60}")
    print(f"[*] Organisms analyzed: {', '.join(sorted(organism_groups))}")
    for label, seqs in sorted(organism_groups.items()):
        run_one_group(label, seqs, OUT_ROOT_ORGANISM, "Organism")

    print(f"\n{'#' * 60}\n# GENUS-LEVEL CpG ISLAND ANALYSIS\n{'#' * 60}")
    print(f"[*] Genera analyzed: {', '.join(sorted(genus_groups))}")
    for label, seqs in sorted(genus_groups.items()):
        run_one_group(label, seqs, OUT_ROOT_GENUS, "Genus")

    print(f"\n[\u2713] Organism-wise + genus-wise CpG island analysis complete.")
    print(f"    Per-organism outputs -> {OUT_ROOT_ORGANISM}/<Organism>/")
    print(f"    Per-genus outputs    -> {OUT_ROOT_GENUS}/<Genus>/")


if __name__ == "__main__":
    main()
