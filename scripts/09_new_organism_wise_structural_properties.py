#!/usr/bin/env python3
"""
STAGE 09-organism — Organism-wise AND Genus-wise DNA Structural Property
Analysis.

WHY GROUPING DOESN'T REUSE 09's OWN extract_clean_species()
---------------------------------------------------------------------------
09's extract_clean_species() only recognizes species via an
'[organism=...]'-style FASTA header tag, or a short hardcoded list of
Plasmodium accession prefixes (PF, PV, PB, PCHAS, PKNH, PY17X, PCOAH). This
pipeline's promoter/gene FASTA headers are plain ">{seq_id}" with no
organism tag, and seq_id is built as "{accession}_{gene}" -- so for any
non-Plasmodium organism (Babesia, Toxoplasma, ...) this heuristic either
falls through to a generic "Plasmodium_sp" bucket, or worse, can
FALSE-POSITIVE match a Plasmodium prefix by coincidence (e.g. a Babesia
gene symbol that happens to start with "PB"). Exactly the same issue
flagged for stage 010 when building the Step 23 comparative matrix.
Grouping here instead goes through family_labels.csv's seq_id -> organism
column, the one join key that's reliable across every stage of this
pipeline regardless of genus.

WHY "Genes" (family_sequences.fasta) NEEDS A LOUD WARNING, NOT A SILENT
NaN WALL
---------------------------------------------------------------------------
family_sequences.fasta is written by stage 01 from UniProt/NCBI PROTEIN
records -- it is amino acid sequence, not DNA. Every structural-property
lookup table in this script (NN_DG, BEND_LOOKUP, MGW_LOOKUP, ROLL_LOOKUP,
PROT_LOOKUP, HELT_LOOKUP) is keyed on DNA dinucleotides. dinuc_profile()
scoring a protein sequence against them doesn't raise an error -- it just
returns NaN for nearly every window, since amino-acid 2-mers essentially
never match a DNA dinucleotide key. That's a pre-existing property of the
original pooled script too, but splitting output into one CSV per
organism/genus would otherwise turn one unexplained oddity into dozens of
silently-all-NaN files with no indication why. This script runs a sniff
test (looks_like_protein()) on every sequence and prints an explicit
warning wherever it's a real risk, rather than leaving you to notice a wall
of NaNs on your own.

BOTH RESOLUTION LEVELS, ONE RUN (matching 08-organism's dual-level pattern)
---------------------------------------------------------------------------
  outputs/DNA_Structural_Properties_By_Organism/<Organism>/<SeqType>_...
  outputs/DNA_Structural_Properties_By_Genus/<Genus>/<SeqType>_...
run for both "Promoters" and "Genes" Seq_Types wherever each file exists.

WHY ORGANISM-LEVEL AND GENUS-LEVEL DON'T GET THE SAME PLOTS
---------------------------------------------------------------------------
The original script's second figure (faceted boxplots by Species) is only
meaningful when a group contains MULTIPLE species to facet across. An
organism-level group is, by definition, one species -- faceting it by
Species would just draw one column. So:
  - ORGANISM level: only the unified all-properties boxplot.
  - GENUS level: the unified boxplot AND the by-species faceted grid
    (scoped to just that genus's own species, which is arguably more
    useful than the original pooled version faceting across every genus
    at once) -- but only when that genus actually has >1 species present;
    otherwise the faceted grid is skipped with a one-line explanation
    rather than silently produced as a trivial single-column plot.

MIN_SEQS_PER_GROUP defaults to 1 (nothing skipped), same reasoning as
08-organism: every row here is an independent, deterministic per-sequence
calculation, not a stochastic control that needs a minimum sample size to
be statistically valid. Raise it yourself if you want to exclude groups
too small to be a meaningful comparison unit.

REUSE
------
The lookup tables (NN_DG, BEND_LOOKUP, MGW_LOOKUP, ROLL_LOOKUP,
PROT_LOOKUP, HELT_LOOKUP), WINDOW, safe_nanmean, dinuc_profile,
curvature_score, nucleosome_score, and summarize_sequence are imported
directly from 09, unchanged. extract_clean_species and locate_sequence_files
are NOT reused, for the reasons above.

Requirements:
    pip install pandas numpy matplotlib seaborn biopython
"""

import os
import sys
import glob
import importlib.util

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE09_FILENAME = None   # set explicitly to skip auto-detection

PROMOTER_FASTA_PATH = "outputs/promoter_sequences.fasta"
GENE_FASTA_PATH = "outputs/family_sequences.fasta"   # protein, per module docstring -- expect warnings
LABELS_PATH = "outputs/family_labels.csv"

OUT_ROOT_ORGANISM = "outputs/DNA_Structural_Properties_By_Organism"
OUT_ROOT_GENUS = "outputs/DNA_Structural_Properties_By_Genus"

MIN_SEQS_PER_GROUP = 1   # see module docstring

PROPERTY_COLS = [
    'Duplex Stability', 'Flexibility', 'Bendability', 'Curvature',
    'Nucleosome Score', 'MGW', 'Roll Angle', 'Propeller Twist', 'Helical Twist'
]
# ===========================================================================


def _find_stage09_module():
    if STAGE09_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE09_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE09_FILENAME is set to {STAGE09_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(
        p for p in glob.glob(os.path.join(_THIS_DIR, "09_*.py"))
        if os.path.basename(p) != os.path.basename(__file__)
    )
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-09 structural properties script "
                  f"(09_*.py) in {_THIS_DIR}. Set STAGE09_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 09_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE09_FILENAME above to override.\n")
    return candidates[0]


_STAGE09_PATH = _find_stage09_module()
_spec = importlib.util.spec_from_file_location("structural_base", _STAGE09_PATH)
sdb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sdb)  # defines functions/constants/lookup tables --
                                 # note: its module-level code ALSO runs
                                 # locate_sequence_files()/the pooled pipeline
                                 # at import time (09 has no __main__ guard),
                                 # so this WILL also produce the original
                                 # pooled outputs as an import side effect.
                                 # That's expected here (unlike 06b/06c/07c/08,
                                 # which do have __main__ guards) -- it's not
                                 # a bug in this script, just how 09 is
                                 # written; the pooled run's outputs are
                                 # harmless to leave in place alongside this
                                 # script's per-group ones.

os.makedirs(OUT_ROOT_ORGANISM, exist_ok=True)
os.makedirs(OUT_ROOT_GENUS, exist_ok=True)


def sanitize_organism_name(name):
    """Same alnum-only sanitizer used throughout the pipeline (01, 05, 06b/
    06c/07c/08-organism, 23_comparative_matrix.py)."""
    import re
    if not name or str(name).strip() == "" or str(name).lower() == "nan":
        return "Unknown_organism"
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_")


def looks_like_protein(seq_str, sample_size=200):
    """Sniff test: DNA/RNA is overwhelmingly A/C/G/T/N/U. A sequence with
    many characters outside that set is almost certainly protein -- flagged
    at >10% non-DNA letters in a sample, well above what ambiguity codes
    alone would ever produce in a real genomic sequence."""
    sample = seq_str[:sample_size].upper()
    if not sample:
        return False
    non_dna = sum(1 for c in sample if c not in "ACGTNU")
    return (non_dna / len(sample)) > 0.10


def load_seq_to_organism():
    if not os.path.exists(LABELS_PATH):
        sys.exit(f"[!] {LABELS_PATH} not found -- run the fetch stage (01) first.")
    labels_df = pd.read_csv(LABELS_PATH)
    if "seq_id" not in labels_df.columns or "organism" not in labels_df.columns:
        sys.exit(f"[!] {LABELS_PATH} is missing 'seq_id' and/or 'organism' columns.")
    return dict(zip(labels_df["seq_id"], labels_df["organism"]))


def build_annotated_records(records, seq_to_org):
    """Returns (annotated, n_unmapped) where annotated is a list of
    (rec, org_label, genus_label) for every record with a resolvable
    organism."""
    annotated = []
    n_unmapped = 0
    for rec in records:
        organism = seq_to_org.get(rec.id)
        if not organism or not str(organism).strip():
            n_unmapped += 1
            continue
        org_label = sanitize_organism_name(organism)
        genus_label = sanitize_organism_name(str(organism).strip().split()[0])
        annotated.append((rec, org_label, genus_label))
    return annotated, n_unmapped


def compute_summary_rows(entries, seq_type_label):
    """entries: list of (rec, org_label) pairs -- org_label becomes that
    row's 'Species' value, which is what makes genus-level faceting-by-
    species meaningful (varies within a genus; constant within an
    organism-level group)."""
    rows = []
    protein_like = 0
    for rec, org_label in entries:
        seq_str = str(rec.seq)
        if not seq_str:
            continue
        if looks_like_protein(seq_str):
            protein_like += 1
        gene_id = rec.id.split()[0]
        rows.append(sdb.summarize_sequence(seq_str, gene_id, org_label, seq_type_label))
    return pd.DataFrame(rows), protein_like


def plot_unified_boxplot(df_melted, group_label, seq_type, out_path):
    """Direct port of 09's own single-canvas boxplot, parameterized by
    group instead of hardcoded to the pooled dataset."""
    plt.figure(figsize=(14, 7), dpi=300)
    sns.boxplot(data=df_melted, x='Property', y='Value', hue='Property',
                palette="Set3", showfliers=False, boxprops=dict(alpha=0.8),
                linewidth=1.2, legend=False)
    sns.stripplot(data=df_melted, x='Property', y='Value',
                  color="black", alpha=0.3, size=3, jitter=0.25)
    plt.title(f"{group_label} — {seq_type} — Comparative DNA Structural Properties Profile",
              fontsize=14, fontweight='bold', pad=15)
    plt.xlabel("Structural Property", fontsize=12, fontweight='bold', labelpad=10)
    plt.ylabel("Value Range / Metric Score", fontsize=12, fontweight='bold', labelpad=10)
    plt.xticks(rotation=25, ha='right', fontsize=11)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_faceted_by_species(df_melted, group_label, seq_type, out_path):
    """Direct port of 09's own faceted-by-species grid, scoped to one
    genus's own species instead of every genus pooled together."""
    g = sns.catplot(data=df_melted, x='Species', y='Value', hue='Species', col='Property',
                     col_wrap=3, kind='box', sharey=False, height=3.5, aspect=1.2,
                     palette='Blues_d', showfliers=False, legend=False)
    g.set_xticklabels(rotation=45, ha='right', fontstyle='italic')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle(f"{group_label} — {seq_type} — All DNA Structural Properties by Species",
                   fontsize=14, fontweight='bold')
    g.savefig(out_path, dpi=300)
    plt.close('all')


def _report_protein_warning(label, seq_type, protein_like, total):
    if protein_like > 0:
        print(f"    [!] {label}/{seq_type}: {protein_like}/{total} sequence(s) look like "
              f"PROTEIN, not DNA -- structural property columns for these rows will be "
              f"mostly/entirely NaN (every lookup table here is keyed on DNA dinucleotides). "
              f"Expected if '{seq_type}' points at family_sequences.fasta (protein, per stage 01) "
              f"rather than a genomic/CDS DNA file.")


def run_organism_group(org_label, entries_by_seqtype):
    group_dir = os.path.join(OUT_ROOT_ORGANISM, org_label)
    for seq_type, entries in entries_by_seqtype.items():
        if len(entries) < MIN_SEQS_PER_GROUP:
            print(f"    [skip] {org_label}/{seq_type}: only {len(entries)} sequence(s) "
                  f"(< MIN_SEQS_PER_GROUP={MIN_SEQS_PER_GROUP}).")
            continue

        df, protein_like = compute_summary_rows(entries, seq_type)
        if df.empty:
            continue
        _report_protein_warning(org_label, seq_type, protein_like, len(entries))

        os.makedirs(group_dir, exist_ok=True)
        csv_path = os.path.join(group_dir, f"{org_label}_{seq_type}_structural_summary.csv")
        df.to_csv(csv_path, index=False)

        df_melted = df.melt(id_vars=['Species', 'Gene_ID'], value_vars=PROPERTY_COLS,
                             var_name='Property', value_name='Value').dropna()
        if df_melted.empty:
            print(f"    [!] {org_label}/{seq_type}: all structural values are NaN -- skipping plot.")
            continue

        plot_path = os.path.join(group_dir, f"{org_label}_{seq_type}_structural_boxplot.png")
        plot_unified_boxplot(df_melted, org_label, seq_type, plot_path)
        print(f"    [\u2713] {org_label}/{seq_type}: {csv_path}")
        print(f"    [\u2713] {org_label}/{seq_type}: {plot_path}")


def run_genus_group(genus_label, entries_by_seqtype):
    group_dir = os.path.join(OUT_ROOT_GENUS, genus_label)
    for seq_type, entries in entries_by_seqtype.items():
        if len(entries) < MIN_SEQS_PER_GROUP:
            print(f"    [skip] {genus_label}/{seq_type}: only {len(entries)} sequence(s) "
                  f"(< MIN_SEQS_PER_GROUP={MIN_SEQS_PER_GROUP}).")
            continue

        df, protein_like = compute_summary_rows(entries, seq_type)
        if df.empty:
            continue
        _report_protein_warning(genus_label, seq_type, protein_like, len(entries))

        os.makedirs(group_dir, exist_ok=True)
        csv_path = os.path.join(group_dir, f"{genus_label}_{seq_type}_structural_summary.csv")
        df.to_csv(csv_path, index=False)

        df_melted = df.melt(id_vars=['Species', 'Gene_ID'], value_vars=PROPERTY_COLS,
                             var_name='Property', value_name='Value').dropna()
        if df_melted.empty:
            print(f"    [!] {genus_label}/{seq_type}: all structural values are NaN -- skipping plots.")
            continue

        plot_path = os.path.join(group_dir, f"{genus_label}_{seq_type}_structural_boxplot.png")
        plot_unified_boxplot(df_melted, genus_label, seq_type, plot_path)
        print(f"    [\u2713] {genus_label}/{seq_type}: {csv_path}")
        print(f"    [\u2713] {genus_label}/{seq_type}: {plot_path}")

        n_species = df['Species'].nunique()
        if n_species > 1:
            facet_path = os.path.join(group_dir, f"{genus_label}_{seq_type}_by_species_grid.png")
            plot_faceted_by_species(df_melted, genus_label, seq_type, facet_path)
            print(f"    [\u2713] {genus_label}/{seq_type}: {facet_path}")
        else:
            print(f"    [i] {genus_label}/{seq_type}: only 1 species present in this genus -- "
                  f"skipping the by-species faceted grid (nothing to facet).")


def main():
    print("=" * 65)
    print("  ORGANISM-WISE + GENUS-WISE DNA STRUCTURAL PROPERTY ANALYSIS  ")
    print("=" * 65)

    seq_to_org = load_seq_to_organism()

    file_targets = []
    if os.path.exists(PROMOTER_FASTA_PATH):
        file_targets.append(("Promoters", PROMOTER_FASTA_PATH))
    if os.path.exists(GENE_FASTA_PATH):
        file_targets.append(("Genes", GENE_FASTA_PATH))
    if not file_targets:
        sys.exit(f"[!] Neither '{PROMOTER_FASTA_PATH}' nor '{GENE_FASTA_PATH}' found.")

    organism_entries, genus_entries = {}, {}

    for seq_type, fpath in file_targets:
        print(f"\n[*] Loading {seq_type} from {fpath} ...")
        records = list(SeqIO.parse(fpath, "fasta"))
        annotated, n_unmapped = build_annotated_records(records, seq_to_org)
        if n_unmapped:
            print(f"    [!] {n_unmapped} sequence(s) had no matching organism in "
                  f"{LABELS_PATH} and were excluded from grouping.")

        protein_like_total = sum(1 for rec, _, _ in annotated if looks_like_protein(str(rec.seq)))
        if protein_like_total > 0:
            print(f"    [!] {seq_type}: {protein_like_total}/{len(annotated)} sequence(s) look "
                  f"like PROTEIN, not DNA. See module docstring -- structural values for these "
                  f"will be mostly/entirely NaN, not an error.")

        for rec, org_label, genus_label in annotated:
            organism_entries.setdefault(org_label, {}).setdefault(seq_type, []).append((rec, org_label))
            genus_entries.setdefault(genus_label, {}).setdefault(seq_type, []).append((rec, org_label))

    print(f"\n{'#' * 60}\n# ORGANISM-LEVEL STRUCTURAL PROPERTY ANALYSIS\n{'#' * 60}")
    print(f"[*] Organisms analyzed: {', '.join(sorted(organism_entries))}")
    for org_label in sorted(organism_entries):
        print(f"\n{'=' * 60}\n[Organism] {org_label}\n{'=' * 60}")
        run_organism_group(org_label, organism_entries[org_label])

    print(f"\n{'#' * 60}\n# GENUS-LEVEL STRUCTURAL PROPERTY ANALYSIS\n{'#' * 60}")
    print(f"[*] Genera analyzed: {', '.join(sorted(genus_entries))}")
    for genus_label in sorted(genus_entries):
        print(f"\n{'=' * 60}\n[Genus] {genus_label}\n{'=' * 60}")
        run_genus_group(genus_label, genus_entries[genus_label])

    print(f"\n[\u2713] Organism-wise + genus-wise structural property analysis complete.")
    print(f"    Per-organism outputs -> {OUT_ROOT_ORGANISM}/<Organism>/")
    print(f"    Per-genus outputs    -> {OUT_ROOT_GENUS}/<Genus>/")


if __name__ == "__main__":
    main()
