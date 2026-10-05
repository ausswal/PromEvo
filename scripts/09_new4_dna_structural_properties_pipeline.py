#!/usr/bin/env python3
"""
DNA STRUCTURAL PROPERTY EVOLUTION PIPELINE (WARNING-FREE PUBLICATION PLOTTING)
Authors: Swarup Das and Subarna Thakur
Affiliation: Department of Bioinformatics, University of North Bengal
"""

import os
import glob
import re
import warnings
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress harmless NumPy/Runtime warnings globally
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ==============================================================================
# 0. CONFIGURATION & FILE LOCATOR (SAVING INTO outputs/)
# ==============================================================================
BASE_OUTPUT = "outputs"
OUT_ROOT = os.path.join(BASE_OUTPUT, "DNA_Structural_Property_Results")
PLOT_DIR = os.path.join(OUT_ROOT, "Plots")

os.makedirs(BASE_OUTPUT, exist_ok=True)
os.makedirs(OUT_ROOT, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

WINDOW = 15  # Sliding window size (bp)


def locate_sequence_files():
    """Locates the PROMOTER FASTA only.

    family_sequences.fasta is intentionally NOT scanned: it holds protein
    sequences, and every lookup table below is a DNA dinucleotide table, so
    protein input produced meaningless (mostly NaN) values. Search order
    matches the rest of the pipeline (outputs/ first); the first hit wins so
    the same file can't be counted twice via '.', '..' and 'outputs'."""
    for r in ['outputs', '.', '..']:
        promoter_f = os.path.join(r, "promoter_sequences.fasta")
        if os.path.exists(promoter_f):
            return [("Promoters", promoter_f)]
    return []


# ==============================================================================
# 1. THERMODYNAMIC & STRUCTURAL LOOKUP TABLES
# ==============================================================================
NN_DG = {
    'AA': -1.00, 'AT': -0.88, 'AC': -1.44, 'AG': -1.28,
    'TA': -0.58, 'TT': -1.00, 'TC': -1.30, 'TG': -1.45,
    'CA': -1.45, 'CT': -1.28, 'CC': -1.84, 'CG': -2.17,
    'GA': -1.30, 'GT': -1.44, 'GC': -2.24, 'GG': -1.84
}

BEND_LOOKUP = {
    'AA': 0.6, 'AT': 0.7, 'AC': 0.4, 'AG': 0.4,
    'TA': 1.0, 'TT': 0.6, 'TC': 0.4, 'TG': 0.4,
    'CA': 0.4, 'CT': 0.4, 'CC': 0.2, 'CG': 0.1,
    'GA': 0.4, 'GT': 0.4, 'GC': 0.1, 'GG': 0.2
}

MGW_LOOKUP = {'AA': 5.4, 'AT': 5.2, 'AC': 5.0, 'AG': 4.8, 'TA': 5.6, 'TT': 5.4, 'TC': 5.1, 'TG': 4.9, 'CA': 4.8, 'CT': 4.6, 'CC': 4.2, 'CG': 3.8, 'GA': 4.7, 'GT': 4.5, 'GC': 3.9, 'GG': 4.1}
ROLL_LOOKUP = {'AA': 1.2, 'AT': 2.0, 'AC': 0.5, 'AG': 1.8, 'TA': 4.5, 'TT': 1.2, 'TC': 1.1, 'TG': 0.5, 'CA': 0.8, 'CT': 1.5, 'CC': -1.2, 'CG': -2.5, 'GA': 1.2, 'GT': 0.8, 'GC': -1.5, 'GG': -0.5}
PROT_LOOKUP = {'AA': -15.2, 'AT': -14.8, 'TC': -12.5, 'TA': -16.0, 'TT': -15.2, 'GC': -7.5, 'CG': -6.8, 'CC': -9.2, 'GG': -10.1, 'AC': -11.0, 'AG': -12.1, 'CA': -11.5, 'CT': -12.0, 'GA': -11.8, 'GT': -11.2, 'TG': -10.8}
HELT_LOOKUP = {'AA': 35.6, 'AT': 33.2, 'TA': 36.0, 'TT': 35.6, 'GC': 31.5, 'CG': 30.2, 'CC': 32.4, 'GG': 32.0, 'AC': 33.8, 'AG': 34.1, 'CA': 33.5, 'CT': 34.0, 'GA': 33.9, 'GT': 34.2, 'TC': 33.7, 'TG': 34.0}

# ==============================================================================
# 2. HELPER FUNCTIONS FOR SAFE MEAN CALCULATION
# ==============================================================================
def safe_nanmean(arr):
    """Calculates mean while ignoring NaNs without throwing empty slice warnings."""
    if arr is None or len(arr) == 0:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        val = np.nanmean(arr)
        return np.nan if np.isnan(val) else val

def dinuc_profile(seq_str, lookup_dict, window=WINDOW):
    seq_str = seq_str.upper()
    if len(seq_str) < 2:
        return np.array([np.nan])
    steps = [seq_str[i:i+2] for i in range(len(seq_str) - 1)]
    vals = np.array([lookup_dict.get(s, np.nan) for s in steps])
    if len(vals) < window or np.isnan(vals).all():
        return vals
    return pd.Series(vals).rolling(window=window, center=True, min_periods=1).mean().values

def curvature_score(seq_str, period=10.5, window=100, step=10):
    n = len(seq_str)
    if n < 30 or window > n:
        return pd.DataFrame(columns=['pos', 'curvature'])
    steps = [seq_str[i:i+2] for i in range(len(seq_str) - 1)]
    atract = np.array([1 if s in ['AA', 'AT', 'TA', 'TT'] else 0 for s in steps])
    if len(atract) < window:
        return pd.DataFrame(columns=['pos', 'curvature'])
        
    lag = int(round(period))
    records = []
    for s_idx in range(0, len(atract) - window + 1, step):
        win = atract[s_idx : s_idx + window]
        if len(win) <= lag:
            continue
        x1, x2 = win[:-lag], win[lag:]
        if np.std(x1) == 0 or np.std(x2) == 0:
            corr = 0.0
        else:
            v = np.corrcoef(x1, x2)[0, 1]
            corr = 0.0 if np.isnan(v) else v
        records.append({'pos': s_idx, 'curvature': corr})
    return pd.DataFrame(records)

def nucleosome_score(seq_str, window=147, step=20):
    n = len(seq_str)
    if n < 30 or window > n:
        return pd.DataFrame(columns=['pos', 'nucleosome_score'])
    chars = list(seq_str.upper())
    at_indicator = np.array([1 if c in ['A', 'T'] else 0 for c in chars])
    lag = 10
    records = []
    
    for s_idx in range(0, n - window + 1, step):
        win_at = at_indicator[s_idx : s_idx + window]
        if len(win_at) <= lag:
            continue
        x1, x2 = win_at[:-lag], win_at[lag:]
        if np.std(x1) == 0 or np.std(x2) == 0:
            periodicity = 0.0
        else:
            v = np.corrcoef(x1, x2)[0, 1]
            periodicity = 0.0 if np.isnan(v) else v
            
        win_chars = chars[s_idx : s_idx + window]
        max_run, current_run = 0, 0
        for c in win_chars:
            if c in ['A', 'T']:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0
        penalty = min(max_run / 20.0, 1.0)
        records.append({'pos': s_idx, 'nucleosome_score': periodicity - penalty})
        
    return pd.DataFrame(records)

def load_organism_lookup():
    """Reads outputs/family_labels.csv's seq_id -> organism mapping -- the
    same authoritative table stage 01 already builds and the rest of this
    pipeline already trusts for species grouping (see the organism-wise
    MSA/tree scripts). Using it here instead of re-guessing from each
    FASTA header's free text is what fixes both problems: headers whose
    organism text is duplicated 2-3x in the description no longer produce
    duplicated species names, and promoter headers that carry no organism
    text at all no longer silently fall back to a wrong hardcoded species."""
    lookup = {}
    for r in ['.', '..', 'outputs']:
        path = os.path.join(r, "family_labels.csv")
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                if 'seq_id' in df.columns and 'organism' in df.columns:
                    for _, row in df.iterrows():
                        sid = str(row['seq_id']).strip()
                        org = str(row['organism']).strip()
                        if sid and org and org.lower() != 'nan':
                            lookup[sid] = org
                    break
            except Exception:
                pass
    return lookup


ORGANISM_LOOKUP = load_organism_lookup()


def _dedupe_repeated_tokens(tokens):
    """Collapses an exact repeated run of tokens -- e.g. ['Toxoplasma',
    'gondii', 'ME49', 'Toxoplasma', 'gondii', 'ME49'] -- down to a single
    copy, regardless of whether the header repeated the organism name 2x,
    3x, or more."""
    n = len(tokens)
    for period in range(1, n // 2 + 1):
        if n % period == 0 and tokens[:period] * (n // period) == tokens:
            return tokens[:period]
    return tokens


def _clean_organism_tokens(raw):
    """Normalizes any raw organism string into a deduplicated token list,
    e.g. 'Babesia bovis T2Bo Babesia bovis T2Bo' -> ['Babesia','bovis','T2Bo'].
    Returns [] if nothing usable is left."""
    raw = re.sub(r'[^a-zA-Z0-9_ ]', '', str(raw)).strip()
    tokens = [t for t in re.split(r'[ _]+', raw) if t]
    return _dedupe_repeated_tokens(tokens)


def _lookup_organism(gene_id):
    """EXACT match on family_labels.csv's seq_id only. The earlier
    suffix-stripping and startswith() prefix fallbacks could silently assign
    the wrong organism to a gene whose ID merely shares a prefix with another
    one; promoter FASTA IDs match family_labels.csv exactly, so no fuzzy
    matching is needed (or wanted). Returns None if the ID is not listed."""
    return ORGANISM_LOOKUP.get(gene_id)


# ID-prefix fallback map, kept for other Apicomplexan projects this
# pipeline may run against (extend as needed). Used only when neither
# family_labels.csv nor any header text yields an organism string.
_ID_PREFIX_SPECIES = [
    (["PF3D7", "PF10", "PF11", "MAL", "PF"], "Plasmodium falciparum"),
    (["PVP01", "PVX", "PV"], "Plasmodium vivax"),
    (["PBANKA", "PB"], "Plasmodium berghei"),
    (["PCHAS", "PCH"], "Plasmodium chabaudi"),
    (["PKNH", "PK"], "Plasmodium knowlesi"),
    (["PY17X", "PY"], "Plasmodium yoelii"),
    (["PCOAH"], "Plasmodium coatneyi"),
]


def resolve_organism_tokens(record_id, record_desc):
    """Single source of truth for 'what organism is this record from', as a
    deduplicated token list (e.g. ['Babesia','bovis','T2Bo']). Everything
    that needs an organism name -- the genus, the genus+species, or the
    full strain-level organism id -- derives from this same token list, in
    this priority order:
      1) family_labels.csv (authoritative -- already used elsewhere in
         this pipeline for organism grouping)
      2) organism=/species= tags in the FASTA header text
      3) a bare "Genus species" pattern anywhere in the header text
      4) a known ID-prefix map (Plasmodium only, currently)
    Returns [] if nothing usable was found anywhere."""
    gene_id = record_id.split()[0]

    organism = _lookup_organism(gene_id)
    if organism:
        tokens = _clean_organism_tokens(organism)
        if tokens:
            return tokens

    full_header = f"{record_id} {record_desc}"
    org_match = re.search(r'\[organism=(.*?)\]', full_header, re.IGNORECASE) or \
                re.search(r'\[species=(.*?)\]', full_header, re.IGNORECASE) or \
                re.search(r'organism=([A-Za-z_][A-Za-z0-9_ ]*)', full_header, re.IGNORECASE)
    if org_match:
        tokens = _clean_organism_tokens(org_match.group(1))
        if tokens:
            return tokens

    bare_match = re.search(r'\b([A-Z][a-z]+)[ _]([a-z]+)(?![a-z])', full_header)
    if bare_match:
        tokens = _clean_organism_tokens(f"{bare_match.group(1)} {bare_match.group(2)}")
        if tokens:
            return tokens

    id_upper = record_id.upper()
    for prefixes, name in _ID_PREFIX_SPECIES:
        if any(id_upper.startswith(p) for p in prefixes):
            return _clean_organism_tokens(name)

    return []


def extract_clean_species(record_id, record_desc):
    """'Genus_species' (strain dropped) -- e.g. 'Toxoplasma_gondii'. Used
    for the combined/top-level plots and for grouping strains together
    within a genus's species-grid."""
    tokens = resolve_organism_tokens(record_id, record_desc)
    if len(tokens) >= 2:
        return f"{tokens[0]}_{tokens[1]}"
    elif len(tokens) == 1:
        return tokens[0]
    # No prior default is safe to guess here -- unlike the old hardcoded
    # "Plasmodium_sp", which silently mislabeled every non-Plasmodium
    # record (that's what produced the all-"Plasmodium_sp" promoter plot).
    return "Unknown_species"


def extract_genus(record_id, record_desc):
    """Just the genus token, e.g. 'Babesia' -- used for the By_Genus output
    tree, which pools every species/strain under one genus together."""
    tokens = resolve_organism_tokens(record_id, record_desc)
    return tokens[0] if tokens else "Unknown_genus"


def extract_full_organism(record_id, record_desc):
    """The full deduplicated organism id, strain included when known --
    e.g. 'Babesia_bovis_T2Bo' -- used for the By_Organism output tree,
    which is the most granular level (one folder per strain/isolate)."""
    tokens = resolve_organism_tokens(record_id, record_desc)
    return "_".join(tokens) if tokens else "Unknown_organism"

def summarize_sequence(seq_str, gene_id, species, genus, organism, label):
    stab = dinuc_profile(seq_str, NN_DG, WINDOW)
    min_dg, max_dg = min(NN_DG.values()), max(NN_DG.values())
    flex_lookup = {k: (v - min_dg) / (max_dg - min_dg) for k, v in NN_DG.items()}
    flex = dinuc_profile(seq_str, flex_lookup, WINDOW)
    
    bend = dinuc_profile(seq_str, BEND_LOOKUP, WINDOW)
    curv = curvature_score(seq_str)
    nucl = nucleosome_score(seq_str)
    
    mgw = dinuc_profile(seq_str, MGW_LOOKUP, WINDOW)
    roll = dinuc_profile(seq_str, ROLL_LOOKUP, WINDOW)
    prot = dinuc_profile(seq_str, PROT_LOOKUP, WINDOW)
    helt = dinuc_profile(seq_str, HELT_LOOKUP, WINDOW)
    
    return {
        'Seq_Type': label,
        'Species': species,
        'Genus': genus,
        'Organism': organism,
        'Gene_ID': gene_id,
        'Duplex Stability': safe_nanmean(stab),
        'Flexibility': safe_nanmean(flex),
        'Bendability': safe_nanmean(bend),
        'Curvature': safe_nanmean(curv['curvature'].values) if not curv.empty else np.nan,
        'Nucleosome Score': safe_nanmean(nucl['nucleosome_score'].values) if not nucl.empty else np.nan,
        'MGW': safe_nanmean(mgw),
        'Roll Angle': safe_nanmean(roll),
        'Propeller Twist': safe_nanmean(prot),
        'Helical Twist': safe_nanmean(helt)
    }

# ==============================================================================
# 3. PIPELINE EXECUTION LOOP
# ==============================================================================
targets = locate_sequence_files()

if not targets:
    print("[ERROR] Could not locate 'promoter_sequences.fasta' (looked in outputs/, ., ..).")
    exit(1)

all_results = []

for label, fpath in targets:
    print(f"\n==== Processing Structural Properties for Target: {label} ({os.path.basename(fpath)}) ====")
    
    try:
        records = list(SeqIO.parse(fpath, "fasta"))
        if not records:
            continue
            
        print(f"  • Found {len(records)} sequence records")
        
        target_summaries = []
        for rec in records:
            seq_str = str(rec.seq)
            if not seq_str:
                continue
            
            species = extract_clean_species(rec.id, rec.description)
            genus = extract_genus(rec.id, rec.description)
            organism = extract_full_organism(rec.id, rec.description)
            res = summarize_sequence(seq_str, rec.id.split()[0], species, genus, organism, label)
            target_summaries.append(res)
            
        if target_summaries:
            df_target = pd.DataFrame(target_summaries)
            target_csv = os.path.join(OUT_ROOT, f"{label}_structural_summary.csv")
            df_target.to_csv(target_csv, index=False)
            all_results.append(df_target)
            
    except Exception as e:
        print(f"  [!] ERROR processing {fpath} - skipped: {e}")

if all_results:
    master_table = pd.concat(all_results, ignore_index=True)
    master_csv = os.path.join(OUT_ROOT, "Master_Structural_Property_Table.csv")
    master_table.to_csv(master_csv, index=False)
    master_table.to_csv("outputs/dna_structural_properties_summary.csv", index=False)

    # ==============================================================================
    # 4. VISUALIZATION
    # ==============================================================================
    property_cols = [
        'Duplex Stability', 'Flexibility', 'Bendability', 'Curvature',
        'Nucleosome Score', 'MGW', 'Roll Angle', 'Propeller Twist', 'Helical Twist'
    ]

    sns.set_theme(style="whitegrid", font="sans-serif")

    def save_unified_boxplot(df_wide, out_paths, title):
        """One boxplot (all properties on the x-axis) + strip overlay,
        saved to every path in out_paths. Returns False (writes nothing) if
        df_wide has no usable rows for the property columns."""
        id_vars = [c for c in ['Species', 'Genus', 'Organism', 'Gene_ID'] if c in df_wide.columns]
        df_melted = df_wide.melt(id_vars=id_vars, value_vars=property_cols,
                                  var_name='Property', value_name='Value').dropna()
        if df_melted.empty:
            return False

        plt.figure(figsize=(14, 7), dpi=300)
        sns.boxplot(data=df_melted, x='Property', y='Value', hue='Property',
                    palette="Set3", showfliers=False, boxprops=dict(alpha=0.8),
                    linewidth=1.2, legend=False)
        sns.stripplot(data=df_melted, x='Property', y='Value',
                      color="black", alpha=0.3, size=3, jitter=0.25)
        plt.title(title, fontsize=14, fontweight='bold', pad=15)
        plt.xlabel("Structural Property", fontsize=12, fontweight='bold', labelpad=10)
        plt.ylabel("Value Range / Metric Score", fontsize=12, fontweight='bold', labelpad=10)
        plt.xticks(rotation=25, ha='right', fontsize=11)
        plt.tight_layout()
        for p in out_paths:
            plt.savefig(p, dpi=300)
        plt.close()
        return True

    def save_facet_grid(df_wide, out_paths, title, group_col='Species'):
        """Per-property facet grid, boxes grouped by group_col, saved to
        every path in out_paths. Returns False if there's nothing to plot."""
        id_vars = [c for c in [group_col, 'Gene_ID'] if c in df_wide.columns]
        df_melted = df_wide.melt(id_vars=id_vars, value_vars=property_cols,
                                  var_name='Property', value_name='Value').dropna()
        if df_melted.empty:
            return False

        g = sns.catplot(
            data=df_melted, x=group_col, y='Value', hue=group_col, col='Property',
            col_wrap=3, kind='box', sharey=False, height=3.5, aspect=1.2,
            palette='Blues_d', showfliers=False, legend=False
        )
        # NOTE: g.set_xticklabels(rotation=45, ...) without explicit label
        # text looks like it should just re-style the existing group names,
        # but FacetGrid grabs each axis's tick-label Text objects before the
        # canvas has ever drawn -- for a categorical plot those are still
        # empty placeholders at that point, so it silently rotates blank
        # labels on every facet. Forcing a draw first, then restyling each
        # axis's real tick labels in place, fixes it.
        g.fig.canvas.draw()
        for ax in g._bottom_axes:  # bottom row only -- avoids overlapping the row above
            ax.tick_params(axis='x', labelbottom=True)
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontstyle='italic')
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle(title, fontsize=14, fontweight='bold')
        for p in out_paths:
            g.savefig(p, dpi=300)
        plt.close()
        return True

    for seq_type, group in master_table.groupby("Seq_Type"):
        save_unified_boxplot(
            group,
            [os.path.join(PLOT_DIR, f"{seq_type}_ALL_PROPERTIES_UNIFIED_BOXPLOT.png"),
             os.path.join(BASE_OUTPUT, f"{seq_type}_ALL_PROPERTIES_UNIFIED_BOXPLOT.png")],
            f"{seq_type} — Comparative DNA Structural Properties Profile"
        )
        save_facet_grid(
            group,
            [os.path.join(PLOT_DIR, f"{seq_type}_ALL_PROPERTIES_FACETED_GRID.png"),
             os.path.join(BASE_OUTPUT, f"{seq_type}_ALL_PROPERTIES_FACETED_GRID.png")],
            f"{seq_type} — All DNA Structural Properties by Species",
            group_col='Species'
        )

    # ==========================================================================
    # 5. BY-GENUS / BY-ORGANISM OUTPUT TREES (additive -- outputs above are
    #    unchanged). One folder per genus (pooling every species/strain of
    #    that genus, plus a species-level facet grid within it), and one
    #    folder per full organism/strain (most granular -- boxplot only,
    #    nothing left to facet by at that level).
    # ==========================================================================
    OUT_BY_GENUS = os.path.join(BASE_OUTPUT, "DNA_Structural_Properties_By_Genus")
    OUT_BY_ORGANISM = os.path.join(BASE_OUTPUT, "DNA_Structural_Properties_By_Organism")

    for genus, genus_group in master_table.groupby("Genus"):
        genus_dir = os.path.join(OUT_BY_GENUS, genus)
        os.makedirs(genus_dir, exist_ok=True)
        for seq_type, sub in genus_group.groupby("Seq_Type"):
            sub.to_csv(os.path.join(genus_dir, f"{genus}_{seq_type}_structural_summary.csv"), index=False)
            save_unified_boxplot(
                sub, [os.path.join(genus_dir, f"{genus}_{seq_type}_structural_boxplot.png")],
                f"{genus} — {seq_type} — Comparative DNA Structural Properties Profile"
            )
            # Grouped by SPECIES: grouping a single genus's data by Genus gave a
            # one-box-per-panel plot. Species is the informative axis here.
            save_facet_grid(
                sub, [os.path.join(genus_dir, f"{genus}_{seq_type}_by_species_grid.png")],
                f"{genus} — {seq_type} — All DNA Structural Properties by Species",
                group_col='Species'
            )

    for organism, org_group in master_table.groupby("Organism"):
        org_dir = os.path.join(OUT_BY_ORGANISM, organism)
        os.makedirs(org_dir, exist_ok=True)
        for seq_type, sub in org_group.groupby("Seq_Type"):
            sub.to_csv(os.path.join(org_dir, f"{organism}_{seq_type}_structural_summary.csv"), index=False)
            save_unified_boxplot(
                sub, [os.path.join(org_dir, f"{organism}_{seq_type}_structural_boxplot.png")],
                f"{organism} — {seq_type} — Comparative DNA Structural Properties Profile"
            )

    print(f"\nBy-genus outputs written to: {OUT_BY_GENUS}/<Genus>/")
    print(f"By-organism outputs written to: {OUT_BY_ORGANISM}/<Organism>/")

print("\n========================================")
print("Pipeline complete!")
print("Outputs folder created at: outputs/DNA_Structural_Property_Results/")
print("Plots generated in: outputs/DNA_Structural_Property_Results/Plots/ and outputs/")
print("========================================\n")