#!/usr/bin/env python3
"""
STAGE 15b — Promoter CpG Island Discovery & Comparative Analysis.

METHOD LABEL FOR YOUR METHODS SECTION: this is the "ADAPTIVE" CpG-island
method (200-bp window, GC cutoff chosen from the dataset's/group's own global
GC%). Stage 08c is the "CLASSIC" method (100-bp window, fixed GC >= 50%,
Obs/Exp >= 0.6). The two can legitimately disagree on the same sequence.

COORDINATES: start/end are promoter-local (0-based, half-open, the form BED
needs). start_tss/end_tss are TSS-relative: each promoter's 3' end sits at
the TSS, so the last base is -1. end_tss is INCLUSIVE (last island base).

Features:
- Global calculation of GC% and AT% prior to sliding-window analysis.
- Adaptive thresholding:
  * Standard Gardiner-Garden parameters (GC >= 50%) if overall GC >= 50%.
  * Adapted AT-rich CpG parameters if overall GC < 50%.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

from pipeline_config import check_promoter_lengths

# ============================== CONFIGURATION ==============================
PROMOTER_FASTA = "outputs/promoter_sequences.fasta"
OUT_DETAILED_CSV = "outputs/cpg_islands_detailed.csv"
OUT_SUMMARY_CSV = "outputs/cpg_island_summary_per_gene.csv"
OUT_BED = "outputs/cpg_islands.bed"
OUT_PLOT = "outputs/cpg_island_distribution.png"

# Standard Parameters (Gardiner-Garden & Frommer)
WINDOW_SIZE = 200
STEP_SIZE = 1
MIN_OBS_EXP = 0.60
# NOTE: because merged islands are built from 200-bp windows, every island is
# automatically >= WINDOW_SIZE bp. MIN_ISLAND_LEN <= WINDOW_SIZE is therefore a
# no-op (the classic Gardiner-Garden minimum of 200 is satisfied by
# construction). Raise it (e.g. 500, Takai-Jones) if you want a real filter.
MIN_ISLAND_LEN = 200
# ===========================================================================


def calculate_metrics(seq_str):
    """Calculates GC%, AT%, and Observed/Expected CpG Ratio for a sequence."""
    seq = seq_str.upper()
    n = len(seq)
    if n == 0:
        return 0.0, 0.0, 0.0

    count_c = seq.count('C')
    count_g = seq.count('G')
    count_a = seq.count('A')
    count_t = seq.count('T')
    count_cg = seq.count('CG')

    # GC and AT percentages
    gc_pct = ((count_c + count_g) / n) * 100.0
    at_pct = ((count_a + count_t) / n) * 100.0

    # Obs/Exp Ratio = (Count(CG) * N) / (Count(C) * Count(G))
    if count_c * count_g == 0:
        obs_exp = 0.0
    else:
        obs_exp = (count_cg * n) / (count_c * count_g)

    return gc_pct, at_pct, obs_exp


def compute_global_nucleotide_content(records):
    """Calculates total GC% and AT% across all records combined."""
    total_len = 0
    total_g = 0
    total_c = 0
    total_a = 0
    total_t = 0

    for rec in records:
        seq = str(rec.seq).upper()
        total_len += len(seq)
        total_g += seq.count('G')
        total_c += seq.count('C')
        total_a += seq.count('A')
        total_t += seq.count('T')

    if total_len == 0:
        return 0.0, 0.0

    global_gc = ((total_g + total_c) / total_len) * 100.0
    global_at = ((total_a + total_t) / total_len) * 100.0
    return global_gc, global_at


def scan_promoter_for_cpg_islands(fasta_path):
    """Scans all promoter sequences in FASTA for CpG islands."""
    if not os.path.exists(fasta_path):
        print(f"[!] Error: Input FASTA file not found at '{fasta_path}'")
        sys.exit(1)

    print(f"[*] Reading promoter sequences from '{fasta_path}'...")
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"[*] Processing {len(records)} promoter sequences...")
    check_promoter_lengths([len(r.seq) for r in records], stage="08")
    if MIN_ISLAND_LEN <= WINDOW_SIZE:
        print(f"[i] MIN_ISLAND_LEN ({MIN_ISLAND_LEN}) <= WINDOW_SIZE ({WINDOW_SIZE}): "
              f"the island-length filter is a no-op (islands are >= {WINDOW_SIZE} bp by construction).")

    # Step 0: Calculate overall GC and AT content across the dataset
    global_gc, global_at = compute_global_nucleotide_content(records)
    print("\n" + "=" * 55)
    print("        GLOBAL NUCLEOTIDE COMPOSITION ANALYSIS        ")
    print("=" * 55)
    print(f" Global GC Content : {global_gc:.2f}%")
    print(f" Global AT Content : {global_at:.2f}%")

    # Select parameters based on overall GC content
    if global_gc >= 50.0:
        min_gc_threshold = 50.0
        mode_label = "Standard Gardiner-Garden Algorithm (GC >= 50%)"
    else:
        # AT-rich adaptation: set threshold relative to baseline or 35%
        min_gc_threshold = max(30.0, round(global_gc + 10.0, 1))
        mode_label = f"AT-Rich Adapted Mode (Lowered GC Cutoff >= {min_gc_threshold}%)"

    print(f" Execution Mode    : {mode_label}")
    print("=" * 55 + "\n")

    all_islands = []
    gene_summaries = []

    for record in records:
        gene_id = record.id
        seq_str = str(record.seq).upper()
        seq_len = len(seq_str)

        overall_gc, overall_at, overall_obsexp = calculate_metrics(seq_str)

        if seq_len < WINDOW_SIZE:
            gene_summaries.append({
                "gene_id": gene_id,
                "promoter_length": seq_len,
                "overall_gc_pct": round(overall_gc, 2),
                "overall_at_pct": round(overall_at, 2),
                "overall_obs_exp": round(overall_obsexp, 3),
                "cpg_island_count": 0,
                "total_cpg_island_bp": 0,
                "has_cpg_island": False
            })
            continue

        # Step 1: Sliding window scan
        raw_hits = []
        for start in range(0, seq_len - WINDOW_SIZE + 1, STEP_SIZE):
            end = start + WINDOW_SIZE
            subseq = seq_str[start:end]
            gc_pct, at_pct, obs_exp = calculate_metrics(subseq)

            if gc_pct >= min_gc_threshold and obs_exp >= MIN_OBS_EXP:
                raw_hits.append((start, end))

        # Step 2: Merge overlapping window hits
        gene_islands = []
        if raw_hits:
            curr_start, curr_end = raw_hits[0]
            for nxt_start, nxt_end in raw_hits[1:]:
                if nxt_start <= curr_end:
                    curr_end = max(curr_end, nxt_end)
                else:
                    island_len = curr_end - curr_start
                    if island_len >= MIN_ISLAND_LEN:
                        gene_islands.append((curr_start, curr_end, island_len))
                    curr_start, curr_end = nxt_start, nxt_end

            island_len = curr_end - curr_start
            if island_len >= MIN_ISLAND_LEN:
                gene_islands.append((curr_start, curr_end, island_len))

        # Step 3: Record detailed island metrics
        total_cpg_bp = 0
        for island_idx, (st, en, ilen) in enumerate(gene_islands, 1):
            island_seq = seq_str[st:en]
            igc, iat, iobsexp = calculate_metrics(island_seq)
            total_cpg_bp += ilen

            all_islands.append({
                "gene_id": gene_id,
                "island_id": f"{gene_id}_CpGI_{island_idx}",
                "start": st,                     # promoter-local, 0-based
                "end": en,                       # promoter-local, half-open
                "start_tss": st - seq_len,       # TSS-relative (last base = -1)
                "end_tss": (en - 1) - seq_len,   # TSS-relative, inclusive
                "length": ilen,
                "gc_pct": round(igc, 2),
                "at_pct": round(iat, 2),
                "obs_exp_ratio": round(iobsexp, 3),
                "cg_count": island_seq.count('CG')
            })

        # Gene level summary
        gene_summaries.append({
            "gene_id": gene_id,
            "promoter_length": seq_len,
            "overall_gc_pct": round(overall_gc, 2),
            "overall_at_pct": round(overall_at, 2),
            "overall_obs_exp": round(overall_obsexp, 3),
            "cpg_island_count": len(gene_islands),
            "total_cpg_island_bp": total_cpg_bp,
            "has_cpg_island": len(gene_islands) > 0
        })

    return pd.DataFrame(all_islands), pd.DataFrame(gene_summaries), min_gc_threshold


def plot_cpg_summary(df_summary, df_islands, min_gc_threshold):
    """Generates a summary plot of CpG islands in promoters."""
    sns.set_theme(style="whitegrid", font="sans-serif")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)

    # Panel A: Overall GC Content vs Obs/Exp Ratio
    sns.scatterplot(
        data=df_summary,
        x="overall_gc_pct",
        y="overall_obs_exp",
        hue="has_cpg_island",
        palette={True: "#D9534F", False: "#4682B4"},
        alpha=0.8,
        s=60,
        ax=ax1
    )
    ax1.axhline(MIN_OBS_EXP, color="red", linestyle="--", linewidth=1, label="Obs/Exp Cutoff (0.6)")
    ax1.axvline(min_gc_threshold, color="orange", linestyle="--", linewidth=1, label=f"GC Cutoff ({min_gc_threshold}%)")
    ax1.set_title("A. Promoter GC Content vs. Observed/Expected CpG Ratio", fontsize=11, fontweight="bold")
    ax1.set_xlabel("Overall GC Content (%)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Observed / Expected CpG Ratio", fontsize=10, fontweight="bold")
    ax1.legend(title="CpG Island Detected", loc="upper left", fontsize=8.5)

    # Panel B: Distribution of Detected CpG Island Lengths
    if not df_islands.empty:
        sns.histplot(
            data=df_islands,
            x="length",
            kde=True,
            color="#2E8B57",
            bins=15,
            edgecolor="black",
            ax=ax2
        )
        ax2.set_title("B. Length Distribution of Identified CpG Islands", fontsize=11, fontweight="bold")
        ax2.set_xlabel("CpG Island Length (bp)", fontsize=10, fontweight="bold")
        ax2.set_ylabel("Frequency", fontsize=10, fontweight="bold")
    else:
        ax2.text(0.5, 0.5, "No CpG Islands Detected", 
                 ha='center', va='center', fontsize=12, fontweight='bold', color='gray')
        ax2.set_title("B. CpG Island Length Distribution", fontsize=11, fontweight="bold")

    plt.suptitle("Promoter CpG Island Discovery — Adaptive Threshold (200-bp window)", fontsize=13, fontweight="bold", y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.savefig(OUT_PLOT, dpi=300)
    plt.close()
    print(f"[✓] Summary plot saved -> {OUT_PLOT}")


def export_bed_file(df_islands, bed_path):
    """Exports detected CpG islands to a standard BED6 file."""
    if df_islands.empty:
        open(bed_path, 'w').close()
        return

    with open(bed_path, 'w') as f:
        for _, row in df_islands.iterrows():
            f.write(f"{row['gene_id']}\t{row['start']}\t{row['end']}\t{row['island_id']}\t{int(row['gc_pct'])}\t+\n")
    print(f"[✓] BED track exported -> {bed_path}")


def main():
    os.makedirs("outputs", exist_ok=True)

    # 1. Perform CpG Island Search
    df_islands, df_summary, min_gc_threshold = scan_promoter_for_cpg_islands(PROMOTER_FASTA)

    # 2. Export Tabular Results
    df_islands.to_csv(OUT_DETAILED_CSV, index=False)
    df_summary.to_csv(OUT_SUMMARY_CSV, index=False)
    export_bed_file(df_islands, OUT_BED)

    print(f"[✓] Detailed CpG Islands saved -> {OUT_DETAILED_CSV}")
    print(f"[✓] Per-gene promoter summary saved -> {OUT_SUMMARY_CSV}")

    # 3. Print Quick Summary Console Stats
    total_genes = len(df_summary)
    genes_with_islands = df_summary['has_cpg_island'].sum() if total_genes > 0 else 0
    print("\n" + "=" * 55)
    print("              CpG ISLAND DISCOVERY SUMMARY              ")
    print("=" * 55)
    print(f" Total Promoter Sequences Analyzed : {total_genes}")
    print(f" Promoters Containing CpG Island(s): {genes_with_islands} ({(genes_with_islands/total_genes*100 if total_genes else 0):.1f}%)")
    print(f" Total CpG Islands Identified      : {len(df_islands)}")
    print("=" * 55 + "\n")

    # 4. Generate Visual Figure
    plot_cpg_summary(df_summary, df_islands, min_gc_threshold)


if __name__ == "__main__":
    main()