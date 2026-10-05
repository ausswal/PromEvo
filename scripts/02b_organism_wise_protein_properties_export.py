#!/usr/bin/env python3
"""
STAGE 08b — Organism- and genus-wise protein physicochemical property
export, built on 02a_protein_properties_for_all.py (imported directly --
no duplication).

WHY THIS IS SAFE TO SPLIT PER-GROUP (unlike 05f's TSS model)
--------------------------------------------------------------
05f could NOT retrain its ML model per genus -- a trained classifier needs
enough independent examples to generalize, and most individual genera
don't have anywhere near enough on their own, so 05f kept training global
and only broke out evaluation by genus. Physicochemical properties are a
different kind of computation: each row is a deterministic calculation
from ONE sequence (molecular weight, pI, aliphatic index, GRAVY, aa
composition, instability index) with no fitting step and no dependence on
any other sequence. Splitting by organism/genus here is just filtering
rows that were already computed once globally -- there's no statistical
reason to keep it global, so unlike 05f, each group's Mean/Std below is a
genuine per-group statistic, not a global model's per-group performance.

Small groups are still flagged rather than silently plotted: below
MIN_N_FOR_PLOT sequences, a group's CSV is still written in full (the raw
numbers are never withheld) but its distribution plot is skipped, since a
histogram/KDE isn't meaningful evidence with only a couple of points.

Requirements:
    pip install biopython pandas numpy matplotlib seaborn
"""

import os
import sys
import glob
import importlib.util

import pandas as pd

# --------------------------------------------------------------------------
# Reuse 02a directly -- no duplication of the property calculations or the
# plot layout (same convention 05f uses to reuse 05b/05c/05d).
# --------------------------------------------------------------------------
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

_STAGE08_PATH = None
_matches = sorted(glob.glob(os.path.join(_THIS_DIR, "02a_*.py")))
if _matches:
    _STAGE08_PATH = _matches[0]
    if len(_matches) > 1:
        print(f"Multiple 02a_*.py scripts found; using {os.path.basename(_matches[0])}.")
if _STAGE08_PATH is None:
    sys.exit(f"Could not find 02a_*.py in {_THIS_DIR} -- this script builds on its "
             f"property calculations and plot layout.")

_spec = importlib.util.spec_from_file_location("protein_props_base", _STAGE08_PATH)
props_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(props_base)

# ============================== CONFIGURATION ==============================
LABELS_PATH = "outputs/family_labels.csv"          # same file 05f reads for genus mapping
ORG_ROOT = "outputs/protein_properties_by_organism"  # matches promoters_by_organism/<Org>/ convention
GENUS_ROOT = "outputs/protein_properties_by_genus"
MIN_N_FOR_PLOT = 5     # below this many sequences, write the CSV but skip the plot
# ===========================================================================


def extract_genus(organism):
    """First word of the organism name -- same convention as
    05_fetch_promoter_sequences_for_all.py / 05f's genus grouping, so
    outputs/*_by_organism/<Genus>/ lines up across every stage."""
    organism = (organism or "").strip()
    if not organism or organism.lower() == "nan":
        return "Unknown_genus"
    return organism.split()[0]


def load_organism_map(labels_path):
    """Returns {seq_id: organism}. Any protein seq_id NOT found in
    family_labels.csv (e.g. a hand-added FASTA entry, or a mismatch with a
    different run) is grouped as 'Unknown_organism' rather than dropped --
    every sequence in family_sequences.fasta still gets a properties row
    somewhere."""
    org_map = {}
    if not os.path.exists(labels_path):
        print(f"[!] {labels_path} not found -- all sequences will be grouped as "
              f"'Unknown_organism' / 'Unknown_genus'.")
        return org_map
    labels_df = pd.read_csv(labels_path)
    for _, row in labels_df.iterrows():
        org_map[str(row["seq_id"])] = str(row.get("organism", "")).strip() or "Unknown_organism"
    return org_map


def sanitize_folder_name(name):
    """Filesystem-safe folder name from an organism/genus string -- same
    convention as the promoter pipeline's sanitize_organism_name."""
    keep = "".join(ch if ch.isalnum() or ch in (" ", "_", "-") else "_" for ch in str(name))
    return "_".join(keep.split())


def export_group(group_df, group_name, group_root, label):
    """Writes one group's CSV, and (if there's enough data) its
    distribution plot, under group_root/<sanitized_group_name>/."""
    group_dir = os.path.join(group_root, sanitize_folder_name(group_name))
    os.makedirs(group_dir, exist_ok=True)

    csv_path = os.path.join(group_dir, "protein_properties.csv")
    group_df.to_csv(csv_path, index=False)

    n = len(group_df)
    if n < MIN_N_FOR_PLOT:
        print(f"  {label:8s} '{group_name}': {n} sequence(s) -> {csv_path} "
              f"(too few for a distribution plot, skipped)")
        return

    plot_path = os.path.join(group_dir, "protein_properties_distributions.png")
    props_base.save_distribution_plot(
        group_df, plot_path,
        f"Physicochemical Property Profiling -- {group_name} (n={n})"
    )
    print(f"  {label:8s} '{group_name}': {n} sequence(s) -> {csv_path}, {plot_path}")


def main():
    os.makedirs("outputs", exist_ok=True)

    seqs = props_base.read_fasta(props_base.INPUT_FASTA)
    if not seqs:
        print(f"Error: Input file '{props_base.INPUT_FASTA}' not found or empty.")
        return

    print(f"Computing protein properties for {len(seqs)} sequences "
          f"({'with' if props_base.HAVE_BIOPYTHON else 'without'} biopython)...")
    df = props_base.compute_properties_df(seqs)
    print(f"[✓] Computed properties for {len(df)}/{len(seqs)} sequences.\n")

    org_map = load_organism_map(LABELS_PATH)
    df["organism"] = df["seq_id"].map(org_map).fillna("Unknown_organism")
    df["genus"] = df["organism"].apply(extract_genus)

    n_unmapped = int((df["organism"] == "Unknown_organism").sum())
    if n_unmapped and org_map:
        print(f"[!] {n_unmapped} sequence(s) in {props_base.INPUT_FASTA} have no matching "
              f"seq_id in {LABELS_PATH} -- grouped as 'Unknown_organism'.\n")

    organism_counts = df["organism"].value_counts().to_dict()
    genus_counts = df["genus"].value_counts().to_dict()
    print(f"Found {len(organism_counts)} organism(s) across {len(genus_counts)} genus/genera: "
          f"{genus_counts}\n")

    print(f"Exporting per-organism summaries under {ORG_ROOT}/<Organism>/ ...")
    for organism, group_df in df.groupby("organism"):
        export_group(group_df.drop(columns=["organism", "genus"]), organism, ORG_ROOT, "Organism")

    print(f"\nExporting per-genus summaries under {GENUS_ROOT}/<Genus>/ ...")
    for genus, group_df in df.groupby("genus"):
        export_group(group_df.drop(columns=["organism", "genus"]), genus, GENUS_ROOT, "Genus")

    # Also refresh the original global CSV + plot (unchanged from 02a's own
    # output), so anything downstream that already reads
    # outputs/protein_properties.csv keeps working untouched.
    df.drop(columns=["organism", "genus"]).to_csv(props_base.OUT_CSV, index=False)
    props_base.save_distribution_plot(
        df, props_base.OUT_PLOT,
        "Physicochemical Property Profiling of Gene Family Proteins"
    )
    print(f"\n[✓] Also refreshed the combined global CSV/plot -> {props_base.OUT_CSV}, "
          f"{props_base.OUT_PLOT}")

    print("\nNow run: python3 09_domain_architecture.py")


if __name__ == "__main__":
    main()
