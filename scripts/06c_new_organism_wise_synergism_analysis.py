#!/usr/bin/env python3
"""
STAGE 06c-organism — Organism-wise (genus-level) Core Promoter Element
Synergism & Co-occurrence Analysis, with per-genus randomized backgrounds
and per-genus outputs.

WHY THE COMPUTATION MUST BE PER-GENUS (same reasoning as 06b-organism)
---------------------------------------------------------------------------
Like 06b, 06c's co-occurrence and spacing statistics are evaluated against a
RANDOMIZED BACKGROUND built by shuffling each real sequence's own
composition. That control is only statistically valid within a
composition-homogeneous group -- GC content varies enormously by genus in
this pipeline, so pooling every organism into one shared shuffle (what 06c
does today) means a pair's real co-occurrence rate in a low-GC genus can be
diluted or masked by a high-GC genus sharing that same pooled background,
producing a misleading synergism signal for BOTH. So each genus gets its
own composition-matched randomized control, computed independently.

OUTPUT: one full report per genus, matching 06b-organism/03b/04b exactly
---------------------------------------------------------------------------
Each genus gets its own complete report + 4-panel figure -- the SAME
content 06c's original run_synergism_analysis() produces (element
multiplicity distribution, pairwise co-occurrence heatmap, TSS-proximal
pair density, inter-element spacing distribution) -- written into
outputs/Core_Promoter_Synergism_By_Organism/<Genus>/. No side-by-side
comparative figure and no single combined report -- if you want to compare
genera afterward, the per-genus report .txt files and the underlying
co-occurrence percentages in each are directly comparable line-by-line.

FIXES PICKED UP FROM 06c (v2): the "MTE-DPE"/"BRE-TATA" pair-name bug (those curves
were never drawn) is fixed via sgb.resolve_pair, and motifs are now scanned on both
strands (as in 06b). DPE/MTE/Inr hits in upstream-only promoters are
sequence-pattern matches judged against the genus-matched shuffle.

REUSE
------
parse_iupac (and the MOTIF_REGEX dict it builds), read_fasta,
shuffle_sequence_composition, scan_motif_positions, and analyze_synergism
are imported directly from 06c, unchanged. The new code is: the
genus-grouping data loader (same one used in 06b-organism), running
shuffle+analyze_synergism PER GENUS instead of once pooled, and per-genus
report/figure output -- the figure and report generation is a straight port
of 06c's own code, just parameterized by genus instead of hardcoded to the
pooled dataset.

Requirements:
    pip install numpy matplotlib seaborn pandas
"""

import os
import sys
import glob
import itertools
import importlib.util

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE06C_FILENAME = None   # set explicitly to skip auto-detection

PROMOTER_FASTA_PATH = "outputs/promoter_sequences.fasta"   # same source 06c already uses
LABELS_PATH = "outputs/family_labels.csv"                   # confirmed "organism" column (01d)

OUT_ROOT = "outputs/Core_Promoter_Synergism_By_Organism"

MIN_SEQS_PER_GENUS = 3   # a composition-matched shuffle control needs at least a
                          # few sequences per genus to mean anything
# ===========================================================================


def _find_stage06c_module():
    if STAGE06C_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE06C_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE06C_FILENAME is set to {STAGE06C_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(
        p for p in glob.glob(os.path.join(_THIS_DIR, "06c_*.py"))
        if os.path.basename(p) != os.path.basename(__file__)
    )
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-06c synergism script (06c_*.py) "
                  f"in {_THIS_DIR}. Set STAGE06C_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 06c_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE06C_FILENAME above to override.\n")
    return candidates[0]


_STAGE06C_PATH = _find_stage06c_module()
_spec = importlib.util.spec_from_file_location("synergism_base", _STAGE06C_PATH)
sgb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sgb)  # only defines functions/constants -- its __main__
                                 # guard keeps its own pooled run from happening
                                 # here as an import side effect.

os.makedirs(OUT_ROOT, exist_ok=True)


# ==============================================================================
# Organism (genus) data loader -- same shape/logic as 06b-organism's, kept
# self-contained here rather than imported since 06c itself has no equivalent
# loader to reuse.
# ==============================================================================
def load_promoters_by_genus():
    """Groups outputs/promoter_sequences.fasta by genus (first word of the
    'organism' column in outputs/family_labels.csv). Returns
    {genus: {seq_id: seq, ...}} -- same dict shape 06c's own functions
    already expect, just one dict per genus instead of one pooled dict."""
    if not os.path.exists(PROMOTER_FASTA_PATH):
        sys.exit(f"[!] {PROMOTER_FASTA_PATH} not found -- run the promoter-fetch stage first.")
    if not os.path.exists(LABELS_PATH):
        sys.exit(f"[!] {LABELS_PATH} not found -- run the fetch stage (01) first "
                 f"(genus grouping needs its 'organism' column).")

    labels_df = pd.read_csv(LABELS_PATH)
    if "organism" not in labels_df.columns:
        sys.exit(f"[!] {LABELS_PATH} has no 'organism' column -- can't group by genus. "
                 f"Actual columns: {list(labels_df.columns)}")
    organism_by_seqid = dict(zip(labels_df["seq_id"], labels_df["organism"]))

    all_seqs = sgb.read_fasta(PROMOTER_FASTA_PATH)
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
# Per-genus report + figure -- a direct port of 06c's own
# run_synergism_analysis(), parameterized by genus instead of hardcoded to
# the pooled dataset. Content and layout are unchanged.
# ==============================================================================
def run_one_genus(genus, seqs):
    n_promoters = len(seqs)
    print(f"\n{'=' * 60}\n[\u27a4] {genus}: {n_promoters} promoter sequence(s)\n{'=' * 60}")

    rand_seqs = sgb.shuffle_sequence_composition(seqs)  # matched to THIS genus's own composition

    print(f"    [*] Analyzing real {genus} promoters...")
    real_synergy = sgb.analyze_synergism(seqs)
    print(f"    [*] Analyzing {genus}-matched randomized background...")
    rand_synergy = sgb.analyze_synergism(rand_seqs)

    genus_dir = os.path.join(OUT_ROOT, genus)
    os.makedirs(genus_dir, exist_ok=True)

    motif_names = list(sgb.CORE_5_MOTIFS.keys())

    # ---------------- 1. Report ----------------
    report_lines = [
        f"{genus.upper()} — CORE PROMOTER ELEMENT SYNERGISM REPORT ({n_promoters} Sequences)",
        "=" * 70,
        "\n1. ELEMENT COMBINATION FREQUENCY PER PROMOTER:",
    ]
    for i in range(6):
        cnt = real_synergy['element_counts'].count(i)
        pct = (cnt / n_promoters) * 100 if n_promoters else 0.0
        report_lines.append(f"  Promoters with {i} core elements: {cnt}/{n_promoters} ({pct:.1f}%)")

    report_lines.append("\n2. PAIRWISE CO-OCCURRENCE PERCENTAGES (REAL vs. RANDOM):")
    for m1, m2 in itertools.combinations(motif_names, 2):
        pair_key = f"{m1}-{m2}"
        r_cnt = real_synergy['matrix'][m1][m2]
        rnd_cnt = rand_synergy['matrix'][m1][m2]
        r_pct = (r_cnt / n_promoters) * 100 if n_promoters else 0.0
        rnd_pct = (rnd_cnt / n_promoters) * 100 if n_promoters else 0.0
        report_lines.append(f"  {pair_key:<12}: Real={r_cnt} ({r_pct:.1f}%) | Random={rnd_cnt} ({rnd_pct:.1f}%)")

    report_path = os.path.join(genus_dir, f"{genus}_promoter_synergism_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
    print("\n".join(f"    {line}" for line in report_lines))
    print(f"    [\u2713] Report -> {report_path}")

    # ---------------- 2. Multi-panel Figure ----------------
    fig = plt.figure(figsize=(16, 11))
    fig.suptitle(f"{genus} — Synergism and Co-occurrence Analysis of Core Promoter Elements",
                 fontsize=15, fontweight='bold')

    # Subplot A: Element Multiplicity Distribution (0 to 5 Elements)
    ax1 = fig.add_subplot(2, 2, 1)
    counts_real = [real_synergy['element_counts'].count(i) for i in range(6)]
    counts_rand = [rand_synergy['element_counts'].count(i) for i in range(6)]
    x = np.arange(6)
    w = 0.35
    ax1.bar(x - w / 2, counts_real, w, label='Real Promoters (Dark)', color='navy')
    ax1.bar(x + w / 2, counts_rand, w, label='Random Background (Light)', color='lightblue')
    ax1.set_xlabel('Number of Core Promoter Elements Present')
    ax1.set_ylabel('Number of Promoters')
    ax1.set_title(f'A. {genus} — Number of Core Elements per Promoter')
    ax1.set_xticks(x)
    ax1.legend()
    ax1.grid(axis='y', linestyle='--', alpha=0.5)

    # Subplot B: Co-occurrence Heatmap Matrix (%)
    ax2 = fig.add_subplot(2, 2, 2)
    heatmap_data = np.zeros((5, 5))
    for i, m1 in enumerate(motif_names):
        for j, m2 in enumerate(motif_names):
            if i != j:
                heatmap_data[i, j] = (real_synergy['matrix'][m1][m2] / n_promoters) * 100 if n_promoters else 0.0
            else:
                heatmap_data[i, j] = np.nan  # Mask diagonal

    sns.heatmap(heatmap_data, annot=True, fmt=".1f", cmap="YlGnBu",
                xticklabels=motif_names, yticklabels=motif_names, ax=ax2,
                cbar_kws={'label': '% Co-occurrence'})
    ax2.set_title(f'B. {genus} — Core Element Pairwise Co-occurrence Matrix (%)')

    # Subplot C: Positional Co-occurrence Density Near TSS
    ax3 = fig.add_subplot(2, 2, 3)
    target_pairs = ['Inr-MTE', 'TATA-Inr', 'TATA-MTE', 'DPE-MTE', 'TATA-BRE', 'TATA-DPE']
    colors = ['green', 'blue', 'purple', 'darkorange', 'red', 'brown']

    any_curve_c = False
    for pair, color in zip(target_pairs, colors):
        key, positions = sgb.resolve_pair(real_synergy['tss_positions'], pair)
        if key and len(positions) > 1 and len(set(positions)) > 1:
            sns.kdeplot(positions, ax=ax3, label=key, color=color, linewidth=2)
            any_curve_c = True
    if not any_curve_c:
        ax3.text(0.5, 0.5, "Not enough pair data\nfor a density curve",
                  ha='center', va='center', transform=ax3.transAxes, fontsize=10, color='gray')

    ax3.axvline(0, color='black', linestyle=':', label='TSS (0)')
    ax3.set_xlim(-1000, 0)
    ax3.set_xlabel('Position Relative to TSS (bp)')
    ax3.set_ylabel('Pair Density Peak')
    ax3.set_title(f'C. {genus} — Pair Co-occurrence Density Peak near TSS')
    ax3.legend(fontsize=8, loc='upper left')
    ax3.grid(True, linestyle=':', alpha=0.6)

    # Subplot D: Pairwise Distance Spacing Distribution (bp)
    ax4 = fig.add_subplot(2, 2, 4)
    any_curve_d = False
    for pair, color in zip(['Inr-MTE', 'TATA-Inr', 'TATA-DPE', 'TATA-BRE'], ['green', 'blue', 'brown', 'red']):
        key, dists = sgb.resolve_pair(real_synergy['distances'], pair)
        if key and len(dists) > 1 and len(set(dists)) > 1:
            sns.kdeplot(dists, ax=ax4, label=f"{key} Spacing", color=color, linewidth=2)
            any_curve_d = True
    if not any_curve_d:
        ax4.text(0.5, 0.5, "Not enough pair data\nfor a density curve",
                  ha='center', va='center', transform=ax4.transAxes, fontsize=10, color='gray')

    ax4.set_xlim(0, 150)
    ax4.set_xlabel('Distance Between Element Pair (bp)')
    ax4.set_ylabel('Density')
    ax4.set_title(f'D. {genus} — Inter-element Distance / Spacing Distribution')
    ax4.legend(fontsize=8)
    ax4.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plot_path = os.path.join(genus_dir, f"{genus}_promoter_synergism_analysis.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"    [\u2713] Figure -> {plot_path}")


def main():
    print("=" * 65)
    print("  ORGANISM-WISE (GENUS) CORE PROMOTER SYNERGISM ANALYSIS  ")
    print("=" * 65)

    genus_seqs = load_promoters_by_genus()
    print(f"[*] Genera analyzed: {', '.join(sorted(genus_seqs))}")
    for g, s in sorted(genus_seqs.items()):
        print(f"      {g}: {len(s)} sequence(s)")

    for genus, seqs in sorted(genus_seqs.items()):
        run_one_genus(genus, seqs)

    print(f"\n[\u2713] Organism-wise core promoter synergism analysis complete. "
          f"One folder per genus under: {OUT_ROOT}/<Genus>/")


if __name__ == "__main__":
    main()
