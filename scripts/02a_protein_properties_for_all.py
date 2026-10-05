#!/usr/bin/env python3
"""
STAGE 08 — Protein Physicochemical Property Analysis & Publication Plotting

Calculates length, molecular weight, theoretical pI, aliphatic index, GRAVY,
amino acid composition, and optional instability index (via BioPython).
Generates publication-quality distribution plots with labeled axes and KDE overlays.
"""

import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

INPUT_FASTA = "outputs/family_sequences.fasta"
OUT_CSV = "outputs/protein_properties.csv"
OUT_PLOT = "outputs/protein_properties_distributions.png"

# Standard average residue masses (Da), one water molecule (18.02) added per peptide
RESIDUE_MASS = {
    "G": 57.05, "A": 71.08, "S": 87.08, "P": 97.12, "V": 99.13, "T": 101.10,
    "C": 103.14, "L": 113.16, "I": 113.16, "N": 114.10, "D": 115.09, "Q": 128.13,
    "K": 128.17, "E": 129.12, "M": 131.19, "H": 137.14, "F": 147.18, "R": 156.19,
    "Y": 163.18, "W": 186.21,
}
WATER = 18.02

KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

PKA = {
    "N_TERM": 9.69, "C_TERM": 2.34, "D": 3.65, "E": 4.25, "C": 8.33,
    "Y": 10.07, "H": 6.00, "K": 10.53, "R": 12.48
}

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    HAVE_BIOPYTHON = True
except ImportError:
    HAVE_BIOPYTHON = False


def read_fasta(path):
    seqs = {}
    name = None
    if not os.path.exists(path):
        return seqs
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            else:
                seqs[name] += line.upper()
    return seqs


def molecular_weight(seq):
    return sum(RESIDUE_MASS.get(aa, 110.0) for aa in seq) + WATER


def gravy(seq):
    return sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq) / len(seq)


def aliphatic_index(seq):
    n = len(seq)
    return 100 * (seq.count("A") + 2.9 * seq.count("V")
                  + 3.9 * (seq.count("I") + seq.count("L"))) / n


def net_charge_at_pH(seq, pH):
    charge = 1 / (1 + 10 ** (pH - PKA["N_TERM"]))
    charge += seq.count("K") / (1 + 10 ** (pH - PKA["K"]))
    charge += seq.count("R") / (1 + 10 ** (pH - PKA["R"]))
    charge += seq.count("H") / (1 + 10 ** (pH - PKA["H"]))
    charge -= 1 / (1 + 10 ** (PKA["C_TERM"] - pH))
    charge -= seq.count("D") / (1 + 10 ** (PKA["D"] - pH))
    charge -= seq.count("E") / (1 + 10 ** (PKA["E"] - pH))
    charge -= seq.count("C") / (1 + 10 ** (PKA["C"] - pH))
    charge -= seq.count("Y") / (1 + 10 ** (PKA["Y"] - pH))
    return charge


def theoretical_pI(seq, tol=0.001):
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        c = net_charge_at_pH(seq, mid)
        if abs(c) < tol:
            return mid
        if c > 0:
            lo = mid
        else:
            hi = mid
    return mid


def aa_composition(seq):
    n = len(seq)
    return {aa: 100 * seq.count(aa) / n for aa in "ACDEFGHIKLMNPQRSTVWY"}


def compute_properties_df(seqs):
    """Given {seq_id: sequence}, computes the same per-sequence
    physicochemical properties main() always has, and returns them as a
    DataFrame. Factored out so 08b (organism/genus-wise export) can reuse
    it verbatim instead of duplicating the calculation."""
    rows = []
    for name, seq in seqs.items():
        clean_seq = "".join(c for c in seq if c in RESIDUE_MASS)
        if len(clean_seq) < 5:
            print(f"  {name}: skipped due to excessive ambiguous residues/short length")
            continue

        row = {
            "seq_id": name,
            "length": len(seq),
            "molecular_weight": round(molecular_weight(clean_seq), 1),
            "theoretical_pI": round(theoretical_pI(clean_seq), 2),
            "aliphatic_index": round(aliphatic_index(clean_seq), 1),
            "gravy": round(gravy(clean_seq), 3),
        }
        if HAVE_BIOPYTHON:
            try:
                row["instability_index"] = round(ProteinAnalysis(clean_seq).instability_index(), 1)
            except Exception:
                row["instability_index"] = None
        row.update({f"aa_pct_{aa}": round(pct, 1) for aa, pct in aa_composition(clean_seq).items()})
        rows.append(row)

    return pd.DataFrame(rows)


# Same 2x2 layout (Molecular Weight / pI / Aliphatic Index / GRAVY) used by
# both the global plot in main() and every per-organism/per-genus plot in
# 08b -- defined once here so the two never drift apart.
PLOT_CONFIGS = [
    {
        "col": "molecular_weight",
        "title": "Molecular Weight Distribution",
        "xlabel": "Molecular Weight (kDa)",
        "color": "#1f77b4",
        "scale": 1e-3  # Convert Da -> kDa for cleaner axis ticks
    },
    {
        "col": "theoretical_pI",
        "title": "Theoretical Isoelectric Point (pI)",
        "xlabel": "Theoretical Isoelectric Point (pH)",
        "color": "#2ca02c",
        "scale": 1.0
    },
    {
        "col": "aliphatic_index",
        "title": "Aliphatic Index Distribution",
        "xlabel": "Aliphatic Index",
        "color": "#ff7f0e",
        "scale": 1.0
    },
    {
        "col": "gravy",
        "title": "GRAVY (Grand Average of Hydropathy)",
        "xlabel": "GRAVY Score",
        "color": "#d62728",
        "scale": 1.0
    }
]


def save_distribution_plot(df, out_plot_path, suptitle):
    """Draws the publication-quality 2x2 distribution figure (MW / pI /
    aliphatic index / GRAVY) for `df` and saves it to out_plot_path.
    Factored out of main() so 08b can render the exact same figure per
    organism/genus, not a simplified copy.

    Handles small-N groups gracefully: KDE is skipped below 5 points (a
    kernel density estimate is not meaningful evidence with a handful of
    values, and can raise on n<2), and a column with zero usable values
    just hides that subplot rather than crashing."""
    sns.set_theme(style="ticks", font_scale=1.0)
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    for ax, cfg in zip(axes.flat, PLOT_CONFIGS):
        data = df[cfg["col"]].dropna() * cfg["scale"] if cfg["col"] in df.columns else pd.Series(dtype=float)
        if data.empty:
            ax.set_visible(False)
            continue

        mean_val = data.mean()
        std_val = data.std() if len(data) > 1 else 0.0
        use_kde = len(data) >= 5
        n_bins = min(22, max(3, len(data)))

        sns.histplot(
            data,
            kde=use_kde,
            ax=ax,
            color=cfg["color"],
            edgecolor="black",
            linewidth=0.8,
            alpha=0.6,
            bins=n_bins
        )

        if cfg["col"] == "gravy":
            ax.axvline(0, color="black", linestyle=":", linewidth=1.2, label="Hydrophobic boundary (0.0)")

        ax.axvline(mean_val, color="red", linestyle="--", linewidth=1.5, label=f"Mean ({mean_val:.2f})")

        ax.set_title(cfg["title"], fontsize=12, fontweight="bold", pad=10)
        ax.set_xlabel(cfg["xlabel"], fontsize=10, fontweight="bold")
        ax.set_ylabel("Number of Proteins (Frequency)", fontsize=10, fontweight="bold")
        ax.grid(axis="y", linestyle="--", alpha=0.5)

        stats_text = f"Mean: {mean_val:.2f}\nStd Dev: {std_val:.2f}\nN = {len(data)}"
        ax.text(
            0.95, 0.92, stats_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85, edgecolor="gray")
        )
        ax.legend(loc="upper left", fontsize=8.5, frameon=True, facecolor="white")

    plt.suptitle(suptitle, fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_plot_path, dpi=300)
    plt.close()


def main():
    os.makedirs("outputs", exist_ok=True)
    seqs = read_fasta(INPUT_FASTA)
    if not seqs:
        print(f"Error: Input file '{INPUT_FASTA}' not found or empty.")
        return

    print(f"Computing protein properties for {len(seqs)} sequences "
          f"({'with' if HAVE_BIOPYTHON else 'without'} biopython)...")

    df = compute_properties_df(seqs)
    df.to_csv(OUT_CSV, index=False)
    print(f"\n[✓] Computed properties for {len(df)}/{len(seqs)} sequences -> {OUT_CSV}")

    save_distribution_plot(df, OUT_PLOT, "Physicochemical Property Profiling of Gene Family Proteins")
    print(f"[✓] Saved publication plot -> {OUT_PLOT}")
    print("\nNow run: python3 09_domain_architecture.py")


if __name__ == "__main__":
    main()