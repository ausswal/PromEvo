#!/usr/bin/env python3
"""
STAGE 06b-organism — Organism-wise (genus-level) Core Promoter Element
Analysis, with per-genus randomized backgrounds and per-genus outputs.

WHY THE COMPUTATION MUST BE PER-GENUS (unlike 06a's grouping-only change)
---------------------------------------------------------------------------
06b's signal-to-noise ratio (SNR) is computed against a RANDOMIZED
BACKGROUND built by shuffling each real sequence's own composition. That
control is only statistically valid within a composition-homogeneous
group. GC content varies enormously by genus in this pipeline (confirmed
empirically in the companion 06a-organism script: ~64% in one genus vs.
~21% in another on real test data). Pooling every organism into ONE shared
randomized background -- what 06b does today -- means a motif's real
enrichment in a low-GC genus could be diluted or masked by a high-GC
genus's contribution to that same pooled shuffle, producing a misleading
SNR for BOTH.

So unlike 06a (where the underlying stats work identically regardless of
grouping, and the "new" part was purely which grouping key fed them), here
the COMPUTATION itself has to change: each genus needs its OWN
composition-matched randomized control.

OUTPUT: one full report per genus, matching 03b/04b exactly
---------------------------------------------------------------------------
Each genus gets its own complete 4-panel figure -- the SAME layout 06b's
original run_plasmodium_promoter_scan() produces (strand distribution,
positional-window distribution, TATA/Inr positional density profiles) --
plus its own summary CSV and individual-hits CSV, all written into
outputs/Core_Promoter_Results_By_Organism/<Genus>/. No side-by-side
comparative figure and no single combined CSV -- if you want to compare
genera afterward, load the per-genus summary CSVs and join on 'motif'.

UPSTREAM-ONLY WINDOWS (v2): positional bins are Distal (-1000..-151), Proximal
(-150..-31) and Core (-30..-1); the old "-30..+20" / "+21..+200" bins could never
be populated by upstream-only promoters. Summary CSV columns are now
distal_hits_1000_to_151 / proximal_hits_150_to_31 / core_hits_30_to_1.
DPE/MTE/Inr stay in the scan; their hits are sequence-pattern matches judged
by SNR against the genus-matched shuffle (see 06b's docstring/figure footnote).

REUSE
------
reverse_complement, parse_iupac (and the MOTIFS dict it builds), read_fasta,
shuffle_sequence_composition, and scan_motifs_and_positions are imported
directly from 06b, unchanged -- no multiprocessing here either, so (like
04b/06a-organism, unlike 03b) nothing needs to be duplicated. The new code
is: the genus-grouping data loader, running scan+shuffle+SNR PER GENUS
instead of once pooled, and per-genus CSV/figure output -- the figure
function itself is a straight port of 06b's own plotting code, just
parameterized by genus instead of hardcoded to the pooled dataset.

Requirements:
    pip install numpy matplotlib pandas
"""

import os
import sys
import csv
import glob
import importlib.util

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE06B_FILENAME = None   # set explicitly to skip auto-detection

PROMOTER_FASTA_PATH = "outputs/promoter_sequences.fasta"   # same source 06b already uses
LABELS_PATH = "outputs/family_labels.csv"                   # confirmed "organism" column (01d)

OUT_ROOT = "outputs/Core_Promoter_Results_By_Organism"

MIN_SEQS_PER_GENUS = 3   # a composition-matched shuffle control needs at least a
                          # few sequences per genus to mean anything
# ===========================================================================


def _find_stage06b_module():
    if STAGE06B_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE06B_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE06B_FILENAME is set to {STAGE06B_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(
        p for p in glob.glob(os.path.join(_THIS_DIR, "06b_*.py"))
        if os.path.basename(p) != os.path.basename(__file__)
    )
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-06b core-promoter script (06b_*.py) "
                  f"in {_THIS_DIR}. Set STAGE06B_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 06b_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE06B_FILENAME above to override.\n")
    return candidates[0]


_STAGE06B_PATH = _find_stage06b_module()
_spec = importlib.util.spec_from_file_location("core_promoter_base", _STAGE06B_PATH)
cpb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cpb)  # only defines functions/constants -- its __main__
                                # guard keeps its own pooled-Plasmodium run from
                                # happening here as an import side effect.

os.makedirs(OUT_ROOT, exist_ok=True)


# ==============================================================================
# NEW: organism (genus) data loader
# ==============================================================================
def load_promoters_by_genus():
    """Groups outputs/promoter_sequences.fasta by genus (first word of the
    'organism' column in outputs/family_labels.csv). Returns
    {genus: {seq_id: seq, ...}} -- same dict shape 06b's own functions
    already expect, just one dict per genus instead of one pooled dict."""
    if not os.path.exists(PROMOTER_FASTA_PATH):
        sys.exit(f"[!] {PROMOTER_FASTA_PATH} not found -- run the promoter-fetch stage first.")
    if not os.path.exists(LABELS_PATH):
        sys.exit(f"[!] {LABELS_PATH} not found -- run the fetch stage (01c/01d) first "
                 f"(genus grouping needs its 'organism' column).")

    labels_df = pd.read_csv(LABELS_PATH)
    if "organism" not in labels_df.columns:
        sys.exit(f"[!] {LABELS_PATH} has no 'organism' column -- can't group by genus. "
                 f"Actual columns: {list(labels_df.columns)}")
    organism_by_seqid = dict(zip(labels_df["seq_id"], labels_df["organism"]))

    all_seqs = cpb.read_fasta(PROMOTER_FASTA_PATH)
    if not all_seqs:
        sys.exit(f"[!] No sequences parsed from {PROMOTER_FASTA_PATH}.")

    genus_seqs = {}
    n_unmapped = 0
    for seq_id, seq in all_seqs.items():
        organism = organism_by_seqid.get(seq_id)
        if not organism or not str(organism).strip():
            n_unmapped += 1
            continue
        genus = str(organism).strip().split()[0]
        genus_seqs.setdefault(genus, {})[seq_id] = seq

    if n_unmapped:
        print(f"[!] {n_unmapped} promoter sequence(s) in {PROMOTER_FASTA_PATH} had no "
              f"matching/non-empty 'organism' in {LABELS_PATH} and were excluded.")

    skipped = {g: len(s) for g, s in genus_seqs.items() if len(s) < MIN_SEQS_PER_GENUS}
    for g, n in skipped.items():
        print(f"[!] Skipping genus '{g}': only {n} sequence(s) (need >= {MIN_SEQS_PER_GENUS} "
              f"for a meaningful composition-matched randomized control).")
        del genus_seqs[g]

    if not genus_seqs:
        sys.exit(f"[!] No genus had >= {MIN_SEQS_PER_GENUS} promoter sequences.")
    return genus_seqs


# ==============================================================================
# Per-genus scan + per-genus randomized background (statistical necessity --
# see module docstring), using 06b's own scan_motifs_and_positions() and
# shuffle_sequence_composition() unchanged.
# ==============================================================================
def run_one_genus(genus, seqs):
    """Scans one genus, writes its CSV pair, and renders its own full 4-panel
    figure -- all into outputs/Core_Promoter_Results_By_Organism/<Genus>/,
    matching 03b/04b's per-organism folder convention exactly."""
    print(f"\n{'=' * 60}\n[\u27a4] {genus}: {len(seqs)} promoter sequence(s)\n{'=' * 60}")

    rand_seqs = cpb.shuffle_sequence_composition(seqs)  # matched to THIS genus's own composition

    print(f"    [*] Scanning real {genus} promoters (+/- strands)...")
    real_res, real_hits = cpb.scan_motifs_and_positions(seqs)
    print(f"    [*] Scanning {genus}-matched randomized background...")
    rand_res, _ = cpb.scan_motifs_and_positions(rand_seqs)

    summary_rows = cpb.build_summary_rows(real_res, rand_res)
    for row in summary_rows:
        print(f"      {row['motif']:<18} real={row['real_total_hits']:<5} rand={row['rand_total_hits']:<5} "
              f"SNR={row['snr']:<6.2f} core(-30..-1)={row[cpb.SUMMARY_COLUMNS[8]]}")

    genus_dir = os.path.join(OUT_ROOT, genus)
    os.makedirs(genus_dir, exist_ok=True)

    summary_path = os.path.join(genus_dir, f"{genus}_core_promoter_summary_metrics.csv")
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"    [\u2713] Summary metrics -> {summary_path}")

    hits_path = os.path.join(genus_dir, f"{genus}_core_promoter_individual_hits.csv")
    hits_cols = ["seq_id", "motif", "strand", "start_pos_0based", "rel_tss_pos_bp", "matched_sequence"]
    hits_df = pd.DataFrame(real_hits)
    hits_df = hits_df[[c for c in hits_cols if c in hits_df.columns]]
    hits_df.to_csv(hits_path, index=False)
    print(f"    [\u2713] Individual motif hits -> {hits_path}")

    plot_path = os.path.join(genus_dir, f"{genus}_core_promoter_analysis.png")
    plot_core_promoter_figure(real_res, len(seqs), genus, plot_path)


def plot_core_promoter_figure(real_res, num_promoters, genus_label, out_path):
    """Thin wrapper: the figure itself now lives in 06b (cpb.plot_core_promoter_figure) so the
    pooled and per-genus figures share ONE implementation and cannot drift apart. Panels use the
    upstream-only windows (Distal -1000..-151 / Proximal -150..-31 / Core -30..-1)."""
    cpb.plot_core_promoter_figure(real_res, num_promoters, genus_label, out_path)
    print(f"    [\u2713] Figure -> {out_path}")


def main():
    print("=" * 65)
    print("  ORGANISM-WISE (GENUS) CORE PROMOTER ELEMENT ANALYSIS  ")
    print("=" * 65)

    genus_seqs = load_promoters_by_genus()
    print(f"[*] Genera analyzed: {', '.join(sorted(genus_seqs))}")
    for g, s in sorted(genus_seqs.items()):
        print(f"      {g}: {len(s)} sequence(s)")

    for genus, seqs in sorted(genus_seqs.items()):
        run_one_genus(genus, seqs)

    print(f"\n[\u2713] Organism-wise core promoter analysis complete. "
          f"One folder per genus under: {OUT_ROOT}/<Genus>/")


if __name__ == "__main__":
    main()
