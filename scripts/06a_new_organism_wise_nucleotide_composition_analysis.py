#!/usr/bin/env python3
"""
STAGE 06a-organism — Organism-wise (genus-level) Nucleotide Frequency &
Windowed Composition Analysis.

DESIGN DIFFERENCE FROM 03b / 04b -- READ BEFORE ASSUMING THIS FRAGMENTS
INTO PER-GENUS FOLDERS THE SAME WAY THOSE DID
-------------------------------------------------------------------------
Composition analysis is inherently COMPARATIVE -- that's what 06a's own
Part A is FOR (pairwise Z-tests comparing GC% between groups). A GC%/AT%
profile computed for one genus in isolation, with nothing to compare it
against, doesn't answer any question this analysis exists to answer. So
"organism-wise" here means:
  - group by GENUS instead of by gene family (matching the same
    genus-level pooling used in 03b/04b), but
  - keep every genus together in ONE set of comparative outputs -- the
    same pairwise significance table, the same faceted/overlay plots,
    just with "Genus" as the grouping variable instead of "Set" (family)
  - NOT fragmented into separate per-genus folders/plots, unlike the tree
    and alignment-visualization stages, where each organism's tree or MSA
    genuinely does stand alone as its own result.

WINDOW COORDINATES: Window_Start_Rel / axis labels now come from 06a's
COVERED_LEN (990 bp = 33 whole 30-nt windows), so the first window starts at
-990, not -1000 (the previous 10 bp mislabelling is fixed in 06a and picked
up here through the import).

REUSE
------
Every computational function -- calculate_gc_at, two_sample_z_test,
benjamini_hochberg, extract_window_features, process_sequence_windows,
process_set_windows -- and every config constant (WINDOW_SIZE,
UPSTREAM_LEN, MONO_BASES, DI_BASES, ...) are imported directly from 06a,
unchanged. 06a has no multiprocessing (like 04b, unlike 03b), so there's
no spawn-pickling risk and nothing needs to be duplicated.

The genuinely NEW code is:
  1. The data loader: 06a's own loader groups sequences by gene-family
     folder name; this one groups by the "organism" column in
     outputs/family_labels.csv instead, collapsed to genus, applied to the
     same outputs/promoter_sequences.fasta 06a's own fallback loader
     already reads.
  2. The pairwise Z-test wiring: 06a's own main() never actually calls
     two_sample_z_test() or benjamini_hochberg() -- both functions are
     fully implemented, but dead code; the docstring promises "pairwise
     two-sample Z-tests" and "Benjamini-Hochberg correction" that the
     script as given never runs. This script is the first place they get
     wired up into real output (a genus x genus GC% significance table).
     Worth fixing in 06a itself too, if you want that comparison for the
     family-level grouping as well -- not done here since it's outside
     what was asked of this script.
  3. Plot titles/labels say "Genus" instead of "Family" where 06a's
     originals hardcoded the latter.

Requirements:
    pip install biopython numpy pandas scipy matplotlib seaborn
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
from Bio import SeqIO

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE06A_FILENAME = None   # set explicitly to skip auto-detection, e.g.
                            # "06a_Nucleotide_Frequency_Windowed_Composition_Analysis.py"

PROMOTER_FASTA_PATH = "outputs/promoter_sequences.fasta"   # same source 06a's own fallback loader uses
LABELS_PATH = "outputs/family_labels.csv"                   # confirmed "organism" column (01d)

OUT_DIR = "outputs/Nucleotide_Frequency_Results_By_Organism"
PLOT_DIR = os.path.join(OUT_DIR, "Plots")

BH_SIGNIFICANCE_ALPHA = 0.05
# ===========================================================================


def _find_stage06a_module():
    if STAGE06A_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE06A_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE06A_FILENAME is set to {STAGE06A_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(
        p for p in glob.glob(os.path.join(_THIS_DIR, "06a_*.py"))
        if os.path.basename(p) != os.path.basename(__file__)
    )
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-06a composition script (06a_*.py) "
                  f"in {_THIS_DIR}. Set STAGE06A_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 06a_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE06A_FILENAME above to override.\n")
    return candidates[0]


_STAGE06A_PATH = _find_stage06a_module()
_spec = importlib.util.spec_from_file_location("comp_base", _STAGE06A_PATH)
comp_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(comp_base)  # only defines functions/constants -- its __main__
                                      # guard keeps its own family-wise run from
                                      # happening here as an import side effect.

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


# ==============================================================================
# NEW: organism (genus) data loader
# ==============================================================================
def load_promoter_data_by_genus():
    """Same promoter FASTA source 06a's own fallback loader already reads,
    grouped by genus (first word of the 'organism' column in
    outputs/family_labels.csv) instead of by gene-family folder. The
    grouping column is still named "Set" so every downstream function
    borrowed from 06a needs zero changes -- it was already grouping-agnostic
    internally, it just happened to always be fed family names before."""
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

    rows = []
    n_unmapped = 0
    for rec in SeqIO.parse(PROMOTER_FASTA_PATH, "fasta"):
        organism = organism_by_seqid.get(rec.id)
        if not organism or not str(organism).strip():
            n_unmapped += 1
            continue
        genus = str(organism).strip().split()[0]
        rows.append({
            "Set": genus,            # deliberately named "Set" -- see docstring above
            "Species": organism,     # kept for reference/labeling
            "Gene_ID": rec.id,
            "Sequence": str(rec.seq),
        })

    if n_unmapped:
        print(f"[!] {n_unmapped} promoter sequence(s) in {PROMOTER_FASTA_PATH} had no "
              f"matching/non-empty 'organism' in {LABELS_PATH} and were excluded.")

    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit(f"[!] No promoter sequences could be mapped to a genus via {LABELS_PATH}.")
    return df


# ==============================================================================
# NEW: pairwise genus x genus Z-tests, wiring up 06a's own (previously
# unused) functions -- see docstring point 2.
# ==============================================================================
def run_pairwise_gc_zscore_tests(per_seq_df, group_col="Set"):
    groups = sorted(per_seq_df[group_col].dropna().unique())
    rows = []
    for g1, g2 in itertools.combinations(groups, 2):
        x1 = per_seq_df.loc[per_seq_df[group_col] == g1, "GC_pct"].to_numpy()
        x2 = per_seq_df.loc[per_seq_df[group_col] == g2, "GC_pct"].to_numpy()
        result = comp_base.two_sample_z_test(x1, x2)
        if result is None:
            print(f"    [!] Skipped {g1} vs {g2}: fewer than 2 usable sequences in one group.")
            continue
        result["Genus1"] = g1
        result["Genus2"] = g2
        rows.append(result)

    if not rows:
        print("    [!] No genus pair had enough sequences for a Z-test.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["p_adj_BH"] = comp_base.benjamini_hochberg(df["p_value"].to_numpy())
    df[f"significant_BH_{BH_SIGNIFICANCE_ALPHA}"] = df["p_adj_BH"] < BH_SIGNIFICANCE_ALPHA
    cols = ["Genus1", "Genus2", "N1", "N2", "Mean1", "Mean2", "Diff", "SE", "Z",
            "p_value", "p_adj_BH", f"significant_BH_{BH_SIGNIFICANCE_ALPHA}"]
    return df[cols]


def main():
    print("=" * 65)
    print("  ORGANISM-WISE (GENUS) NUCLEOTIDE FREQUENCY & COMPOSITION  ")
    print("=" * 65)

    df_promoters = load_promoter_data_by_genus()
    genera = sorted(df_promoters["Set"].unique())
    print(f"[*] Analyzing promoter sequences for genera: {', '.join(genera)}")
    for g in genera:
        print(f"      {g}: {(df_promoters['Set'] == g).sum()} sequence(s)")

    # --------------------------------------------------------------------------
    # PART A: Whole-sequence GC%/AT%, per sequence -- reuses calculate_gc_at()
    # unchanged from 06a.
    # --------------------------------------------------------------------------
    stats_list = []
    for _, row in df_promoters.iterrows():
        gc_pct, at_pct, length = comp_base.calculate_gc_at(row["Sequence"])
        stats_list.append({
            "Genus": row["Set"], "Species": row["Species"], "Gene_ID": row["Gene_ID"],
            "Length": length, "GC_pct": gc_pct, "AT_pct": at_pct,
        })
    per_seq_df = pd.DataFrame(stats_list)
    per_seq_df.to_csv(os.path.join(OUT_DIR, "Per_Sequence_GC_AT_Content.csv"), index=False)

    print("\n[*] Running pairwise genus x genus GC% Z-tests "
          "(06a defines this but never calls it -- wired up here)...")
    zscore_df = run_pairwise_gc_zscore_tests(per_seq_df.rename(columns={"Genus": "Set"}))
    if not zscore_df.empty:
        zscore_path = os.path.join(OUT_DIR, "Pairwise_Genus_GC_ZTests_BH_corrected.csv")
        zscore_df.to_csv(zscore_path, index=False)
        n_sig = int(zscore_df[f"significant_BH_{BH_SIGNIFICANCE_ALPHA}"].sum())
        print(f"    [\u2713] {len(zscore_df)} genus pair(s) tested, {n_sig} significant "
              f"at BH-adjusted p < {BH_SIGNIFICANCE_ALPHA} -> {zscore_path}")

    # 1. Whole-Sequence Combined AT vs GC Boxplot, by genus
    plt.figure(figsize=(9, 6))
    melted_seq = pd.melt(per_seq_df, id_vars=["Genus"], value_vars=["AT_pct", "GC_pct"],
                          var_name="Metric", value_name="Percentage")
    melted_seq["Metric"] = melted_seq["Metric"].map({"AT_pct": "AT %", "GC_pct": "GC %"})
    sns.boxplot(data=melted_seq, x="Genus", y="Percentage", hue="Metric",
                palette={"AT %": "#d95f02", "GC %": "#7570b3"})
    plt.title("Promoter Whole-Sequence AT% vs GC% Content — Organism (Genus) Comparison")
    plt.xlabel("")
    plt.ylabel("Percentage (%)")
    plt.ylim(0, 100)
    plt.xticks(rotation=35, ha="right")
    plt.legend(title="Composition", loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "Whole_Sequence_AT_vs_GC_By_Genus.png"), dpi=150)
    plt.close()

    # --------------------------------------------------------------------------
    # PART B: Sliding-window composition, by genus -- reuses process_set_windows()
    # unchanged from 06a.
    # --------------------------------------------------------------------------
    window_tables = []
    for genus in genera:
        seqs = df_promoters.loc[df_promoters["Set"] == genus, "Sequence"].tolist()
        wt = comp_base.process_set_windows(seqs)
        wt["Genus"] = genus
        window_tables.append(wt)

    window_df = pd.concat(window_tables, ignore_index=True)
    window_df["AT_pct"] = (window_df["AT_count"] / comp_base.WINDOW_SIZE) * 100.0
    window_df["GC_pct"] = (window_df["GC_count"] / comp_base.WINDOW_SIZE) * 100.0
    window_df.to_csv(os.path.join(OUT_DIR, "Windowed_MonoDi_Nucleotide_Composition.csv"), index=False)

    # 2. Windowed AT% vs GC% Overlay Profiles, faceted by genus
    window_long = pd.melt(window_df, id_vars=["Genus", "Window_Start_Rel"],
                           value_vars=["AT_pct", "GC_pct"],
                           var_name="Nucleotide_Type", value_name="Percentage")
    window_long["Nucleotide_Type"] = window_long["Nucleotide_Type"].map(
        {"AT_pct": "AT %", "GC_pct": "GC %"})

    g = sns.FacetGrid(window_long, col="Genus",
                       hue="Nucleotide_Type", palette={"AT %": "#e66101", "GC %": "#5e3c99"},
                       col_wrap=3, height=4, sharey=True)
    g.map(sns.lineplot, "Window_Start_Rel", "Percentage", marker="o", linewidth=2)
    g.add_legend(title="Metric")
    g.set_axis_labels(f"Position relative to gene start (-{comp_base.COVERED_LEN} to -1)",
                       "Mean Percentage (%)")
    g.set_titles(col_template="Genus: {col_name}")
    g.fig.subplots_adjust(top=0.88)
    g.fig.suptitle(f"Windowed AT% vs GC% Profile — Organism (Genus) Comparison "
                    f"({comp_base.WINDOW_SIZE}-nt Windows)", fontsize=14)
    g.savefig(os.path.join(PLOT_DIR, "Windowed_AT_vs_GC_Overlay_Profile_By_Genus.png"), dpi=150)
    plt.close()

    # 3. Stacked comparative profile, one panel per genus, all in one figure
    fig, axes = plt.subplots(len(genera), 1, figsize=(10, 3 * len(genera)), sharex=True, squeeze=False)
    for idx, genus in enumerate(genera):
        gdata = window_df[window_df["Genus"] == genus].sort_values("Window_Start_Rel")
        ax = axes[idx, 0]
        ax.plot(gdata["Window_Start_Rel"], gdata["AT_pct"], label="AT %", color="#d95f02",
                linewidth=2, marker="o")
        ax.plot(gdata["Window_Start_Rel"], gdata["GC_pct"], label="GC %", color="#7570b3",
                linewidth=2, marker="s")
        ax.axhline(50, color="gray", linestyle="--", alpha=0.5)
        ax.set_ylabel("Percentage (%)")
        ax.set_title(f"Genus: {genus}")
        ax.set_ylim(0, 100)
        ax.legend(loc="upper right")
        ax.grid(True, linestyle=":", alpha=0.6)

    axes[-1, 0].set_xlabel(f"Position relative to anchor (-{comp_base.COVERED_LEN} to -1)")
    plt.suptitle("Comparative AT% and GC% Windowed Profiles Across Organisms (Genus)",
                 y=1.01, fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "Windowed_AT_vs_GC_Comparative_Stack_By_Genus.png"), dpi=150)
    plt.close()

    print("\n========================================")
    print(f"[\u2713] Organism-wise comparative analysis complete. Results saved in: {OUT_DIR}/")
    print("  \u2022 Per_Sequence_GC_AT_Content.csv")
    print("  \u2022 Pairwise_Genus_GC_ZTests_BH_corrected.csv")
    print("  \u2022 Windowed_MonoDi_Nucleotide_Composition.csv")
    print("  \u2022 Plots/Whole_Sequence_AT_vs_GC_By_Genus.png")
    print("  \u2022 Plots/Windowed_AT_vs_GC_Overlay_Profile_By_Genus.png")
    print("  \u2022 Plots/Windowed_AT_vs_GC_Comparative_Stack_By_Genus.png")
    print("========================================\n")


if __name__ == "__main__":
    main()
