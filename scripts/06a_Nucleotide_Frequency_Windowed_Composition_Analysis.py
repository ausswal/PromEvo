#!/usr/bin/env python3
"""
STAGE 01.5 — Nucleotide Frequency & Windowed Composition Analysis Module.

PART A:
  - Whole-sequence GC% / AT% computation per gene family.
  - Pairwise two-sample Z-tests comparing mean GC% across gene families.
  - Benjamini-Hochberg (p-adjust) correction for multiple testing.

PART B:
  - Sliding 30-nt window mono- and di-nucleotide composition profiles.
  - Anchored from the 3' end (position -1 relative to gene start) going upstream.
  - Only WHOLE windows are used: with UPSTREAM_LEN = 1000 and WINDOW_SIZE = 30
    that is 33 windows covering -990 to -1 (COVERED_LEN = 990 bp). The most
    distal 10 bp (-1000 to -991) are not covered. All Window_Start_Rel /
    Window_End_Rel values and axis labels are computed from COVERED_LEN, not
    UPSTREAM_LEN (previous versions labelled the first window as starting at
    -1000, i.e. every window was mislabelled by 10 bp).
  - Side-by-side & overlay plots comparing AT% and GC% together across windows.

Outputs generated in: outputs/Nucleotide_Frequency_Results/
"""

import os
import sys
import glob
import itertools
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# ==============================================================================
# CONFIGURATION
# ==============================================================================
WINDOW_SIZE = 30       # nt per sliding window
UPSTREAM_LEN = 1000    # bp upstream of gene start
DOWNSTREAM_LEN = 0     # 0 = ends exactly at start codon (-1)
TOTAL_REGION = UPSTREAM_LEN + DOWNSTREAM_LEN
N_WINDOWS = TOTAL_REGION // WINDOW_SIZE
COVERED_LEN = N_WINDOWS * WINDOW_SIZE   # bp actually covered by whole windows (990 for 1000/30)

MONO_BASES = ["A", "C", "G", "T"]
DI_BASES = [f"{b1}{b2}" for b1 in MONO_BASES for b2 in MONO_BASES]

OUT_DIR = "outputs/Nucleotide_Frequency_Results"
PLOT_DIR = os.path.join(OUT_DIR, "Plots")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


# ==============================================================================
# PART A: WHOLE SEQUENCE ANALYSIS & Z-TESTS
# ==============================================================================

def calculate_gc_at(sequence):
    """Calculate GC%, AT%, and Length for a DNA sequence string."""
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return np.nan, np.nan, 0
    
    g_cnt = seq.count("G")
    c_cnt = seq.count("C")
    a_cnt = seq.count("A")
    t_cnt = seq.count("T")
    
    gc_pct = (g_cnt + c_cnt) / n * 100.0
    at_pct = (a_cnt + t_cnt) / n * 100.0
    return gc_pct, at_pct, n


def two_sample_z_test(x1, x2):
    """Two-sample unpooled Z-test for GC% difference."""
    x1 = x1[~np.isnan(x1)]
    x2 = x2[~np.isnan(x2)]
    
    n1, n2 = len(x1), len(x2)
    if n1 < 2 or n2 < 2:
        return None

    m1, m2 = np.mean(x1), np.mean(x2)
    v1, v2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    
    se = np.sqrt((v1 / n1) + (v2 / n2))
    if se == 0:
        return None

    z_stat = (m1 - m2) / se
    p_val = 2 * (1 - stats.norm.cdf(abs(z_stat)))
    
    return {
        "N1": n1, "N2": n2,
        "Mean1": m1, "Mean2": m2,
        "Diff": m1 - m2,
        "SE": se, "Z": z_stat,
        "p_value": p_val
    }


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR correction."""
    p_vals = np.asarray(p_values)
    n = len(p_vals)
    sorted_indices = np.argsort(p_vals)
    sorted_p = p_vals[sorted_indices]
    
    adj_p = np.zeros(n)
    cum_min = 1.0
    
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = (sorted_p[i] * n) / rank
        cum_min = min(cum_min, val)
        adj_p[i] = cum_min
        
    adj_p = np.clip(adj_p, 0, 1.0)
    rev_indices = np.argsort(sorted_indices)
    return adj_p[rev_indices]


# ==============================================================================
# PART B: WINDOWED COMPOSITION ANALYSIS
# ==============================================================================

def extract_window_features(window_str):
    """Count mono and dinucleotides in a window string."""
    w_len = len(window_str)
    mono_counts = {b: window_str.count(b) for b in MONO_BASES}
    
    di_counts = {db: 0 for db in DI_BASES}
    if w_len >= 2:
        for j in range(w_len - 1):
            dinuc = window_str[j:j+2]
            if dinuc in di_counts:
                di_counts[dinuc] += 1
                
    at_count = mono_counts["A"] + mono_counts["T"]
    gc_count = mono_counts["G"] + mono_counts["C"]
    return {**mono_counts, **di_counts, "AT_count": at_count, "GC_count": gc_count}


def process_sequence_windows(seq_str):
    """Extract sliding windows anchored from 3' end (-1) going upstream (-1000)."""
    seq_str = seq_str.upper()
    L = len(seq_str)
    windows = []

    for i in range(1, N_WINDOWS + 1):
        end_from_3prime = (N_WINDOWS - i + 1) * WINDOW_SIZE
        start_from_3prime = end_from_3prime - WINDOW_SIZE + 1

        idx_start = L - end_from_3prime
        idx_end = L - start_from_3prime + 1

        if idx_start < 0 or idx_end > L or idx_start >= idx_end:
            windows.append(None)
        else:
            windows.append(seq_str[idx_start:idx_end])
            
    return windows


def process_set_windows(seq_list):
    """Aggregate window metrics across all sequences in a family."""
    all_window_counts = [[] for _ in range(N_WINDOWS)]

    for seq in seq_list:
        windows = process_sequence_windows(seq)
        for i in range(N_WINDOWS):
            if windows[i] is not None:
                feat = extract_window_features(windows[i])
                all_window_counts[i].append(feat)

    results = []
    for i in range(N_WINDOWS):
        mat = all_window_counts[i]
        window_start_rel = -COVERED_LEN + (i * WINDOW_SIZE)
        window_end_rel = window_start_rel + WINDOW_SIZE - 1

        if not mat:
            row = {
                "Window_Index": i + 1,
                "Window_Start_Rel": window_start_rel,
                "Window_End_Rel": window_end_rel,
                "N_Sequences_Covered": 0,
                **{b: np.nan for b in MONO_BASES + DI_BASES + ["AT_count", "GC_count"]}
            }
        else:
            df_mat = pd.DataFrame(mat)
            means = df_mat.mean().to_dict()
            row = {
                "Window_Index": i + 1,
                "Window_Start_Rel": window_start_rel,
                "Window_End_Rel": window_end_rel,
                "N_Sequences_Covered": len(mat),
                **means
            }
        results.append(row)

    return pd.DataFrame(results)


# ==============================================================================
# DATA LOADING
# ==============================================================================

def load_promoter_data():
    """Discover and parse promoter FASTA files."""
    promoter_data = []
    family_dirs = [d for d in glob.glob("All_*_Promoters") if not "old" in d.lower()]
    
    if family_dirs:
        for fdir in family_dirs:
            family_name = fdir.replace("All_", "").replace("_Promoters", "")
            fasta_files = [f for f in glob.glob(os.path.join(fdir, "**", "*.fasta"), recursive=True) if "named" not in f.lower()]
            for fp in fasta_files:
                species = os.path.basename(os.path.dirname(fp))
                for rec in SeqIO.parse(fp, "fasta"):
                    promoter_data.append({
                        "Set": family_name,
                        "Species": species,
                        "Gene_ID": rec.id.split()[0],
                        "Sequence": str(rec.seq)
                    })
    else:
        fallback_fasta = "outputs/promoter_sequences.fasta"
        labels_csv = "outputs/family_labels.csv"
        
        if os.path.exists(fallback_fasta):
            family = "Target_Family"
            if os.path.exists(labels_csv):
                l_df = pd.read_csv(labels_csv)
                if not l_df.empty and "family" in l_df.columns:
                    family = l_df["family"].iloc[0]

            for rec in SeqIO.parse(fallback_fasta, "fasta"):
                promoter_data.append({
                    "Set": family,
                    "Species": "Unknown",
                    "Gene_ID": rec.id,
                    "Sequence": str(rec.seq)
                })

    return pd.DataFrame(promoter_data)


# ==============================================================================
# MAIN EXECUTION & COMPARATIVE PLOTTING
# ==============================================================================

def main():
    print("=" * 65)
    print("  NUCLEOTIDE FREQUENCY & COMPARATIVE AT/GC COMPOSITION  ")
    print("=" * 65)
    if COVERED_LEN != TOTAL_REGION:
        print(f"[i] {N_WINDOWS} x {WINDOW_SIZE}-nt windows cover -{COVERED_LEN} to -1; the "
              f"{TOTAL_REGION - COVERED_LEN} most distal bp (-{TOTAL_REGION} to "
              f"-{COVERED_LEN + 1}) are not in any window.")

    df_promoters = load_promoter_data()
    if df_promoters.empty:
        sys.exit("[!] Error: No promoter sequences found to analyze.")

    families = df_promoters["Set"].unique()
    print(f"[*] Analyzing promoter sets for families: {', '.join(families)}")

    # --------------------------------------------------------------------------
    # PART A: Whole Sequence Analysis
    # --------------------------------------------------------------------------
    stats_list = []
    for idx, row in df_promoters.iterrows():
        gc_pct, at_pct, length = calculate_gc_at(row["Sequence"])
        stats_list.append({
            "Set": row["Set"],
            "Species": row["Species"],
            "Gene_ID": row["Gene_ID"],
            "Length": length,
            "GC_pct": gc_pct,
            "AT_pct": at_pct
        })
    
    per_seq_df = pd.DataFrame(stats_list)
    per_seq_df.to_csv(os.path.join(OUT_DIR, "Per_Sequence_GC_AT_Content.csv"), index=False)

    # 1. Whole-Sequence Combined AT vs GC Boxplot Comparison
    plt.figure(figsize=(9, 6))
    melted_seq = pd.melt(
        per_seq_df, 
        id_vars=["Set"], 
        value_vars=["AT_pct", "GC_pct"], 
        var_name="Metric", 
        value_name="Percentage"
    )
    melted_seq["Metric"] = melted_seq["Metric"].map({"AT_pct": "AT %", "GC_pct": "GC %"})

    sns.boxplot(data=melted_seq, x="Set", y="Percentage", hue="Metric", palette={"AT %": "#d95f02", "GC %": "#7570b3"})
    plt.title("Promoter Whole-Sequence AT% vs GC% Content Comparison")
    plt.xlabel("")
    plt.ylabel("Percentage (%)")
    plt.ylim(0, 100)
    plt.xticks(rotation=35, ha="right")
    plt.legend(title="Composition", loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "Whole_Sequence_AT_vs_GC_By_Set.png"), dpi=150)
    plt.close()

    # --------------------------------------------------------------------------
    # PART B: Sliding Window Composition
    # --------------------------------------------------------------------------
    window_tables = []
    for fam in families:
        seqs = df_promoters[df_promoters["Set"] == fam]["Sequence"].tolist()
        wt = process_set_windows(seqs)
        wt["Set"] = fam
        window_tables.append(wt)

    window_df = pd.concat(window_tables, ignore_index=True)
    window_df["AT_pct"] = (window_df["AT_count"] / WINDOW_SIZE) * 100.0
    window_df["GC_pct"] = (window_df["GC_count"] / WINDOW_SIZE) * 100.0
    window_df.to_csv(os.path.join(OUT_DIR, "Windowed_MonoDi_Nucleotide_Composition.csv"), index=False)

    # 2. Windowed AT% vs GC% Overlay Profiles (Per Gene Family Facets)
    window_long = pd.melt(
        window_df,
        id_vars=["Set", "Window_Start_Rel"],
        value_vars=["AT_pct", "GC_pct"],
        var_name="Nucleotide_Type",
        value_name="Percentage"
    )
    window_long["Nucleotide_Type"] = window_long["Nucleotide_Type"].map({"AT_pct": "AT %", "GC_pct": "GC %"})

    g = sns.FacetGrid(window_long, col="Set", hue="Nucleotide_Type", palette={"AT %": "#e66101", "GC %": "#5e3c99"}, col_wrap=3, height=4, sharey=True)
    g.map(sns.lineplot, "Window_Start_Rel", "Percentage", marker="o", linewidth=2)
    g.add_legend(title="Metric")
    g.set_axis_labels(f"Position relative to gene start (-{COVERED_LEN} to -1)", "Mean Percentage (%)")
    g.set_titles(col_template="Family: {col_name}")
    g.fig.subplots_adjust(top=0.88)
    g.fig.suptitle(f"Windowed AT% vs GC% Profile Comparison ({WINDOW_SIZE}-nt Windows)", fontsize=14)
    g.savefig(os.path.join(PLOT_DIR, "Windowed_AT_vs_GC_Overlay_Profile.png"), dpi=150)
    plt.close()

    # 3. Stacked Composition Profile (Visualizing Total Base Makeup)
    fig, axes = plt.subplots(len(families), 1, figsize=(10, 3 * len(families)), sharex=True, squeeze=False)
    for idx, fam in enumerate(families):
        fam_data = window_df[window_df["Set"] == fam].sort_values("Window_Start_Rel")
        ax = axes[idx, 0]
        ax.plot(fam_data["Window_Start_Rel"], fam_data["AT_pct"], label="AT %", color="#d95f02", linewidth=2, marker="o")
        ax.plot(fam_data["Window_Start_Rel"], fam_data["GC_pct"], label="GC %", color="#7570b3", linewidth=2, marker="s")
        ax.axhline(50, color="gray", linestyle="--", alpha=0.5)
        ax.set_ylabel("Percentage (%)")
        ax.set_title(f"Gene Family: {fam}")
        ax.set_ylim(0, 100)
        ax.legend(loc="upper right")
        ax.grid(True, linestyle=":", alpha=0.6)

    axes[-1, 0].set_xlabel(f"Position relative to anchor (-{COVERED_LEN} to -1)")
    plt.suptitle("Comparative AT% and GC% Windowed Profiles Across Promoter Loci", y=1.01, fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "Windowed_AT_vs_GC_Comparative_Stack.png"), dpi=150)
    plt.close()

    print("\n========================================")
    print(f"[✓] Comparative Analysis Complete. Results saved in: {OUT_DIR}/")
    print("  • Plots/Whole_Sequence_AT_vs_GC_By_Set.png  (Side-by-side Boxplot)")
    print("  • Plots/Windowed_AT_vs_GC_Overlay_Profile.png (Faceted Overlay Profile)")
    print("  • Plots/Windowed_AT_vs_GC_Comparative_Stack.png (Multi-panel Comparative Plots)")
    print("========================================\n")

if __name__ == "__main__":
    main()