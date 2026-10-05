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
    """Locates target FASTA sequence files in current dir, parent dir, or 'outputs/'."""
    targets = []
    possible_roots = ['.', '..', 'outputs']

    for r in possible_roots:
        promoter_f = os.path.join(r, "promoter_sequences.fasta")
        gene_f = os.path.join(r, "family_sequences.fasta")
        
        if os.path.exists(promoter_f) and ("Promoters", promoter_f) not in targets:
            targets.append(("Promoters", promoter_f))
        if os.path.exists(gene_f) and ("Genes", gene_f) not in targets:
            targets.append(("Genes", gene_f))

    return targets


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

def extract_clean_species(record_id, record_desc):
    full_header = f"{record_id} {record_desc}"
    org_match = re.search(r'\[organism=(.*?)\]', full_header, re.IGNORECASE) or \
                re.search(r'\[species=(.*?)\]', full_header, re.IGNORECASE) or \
                re.search(r'organism=(.*?)(?:\s|\]|$)', full_header, re.IGNORECASE)
    
    if org_match:
        sp_name = org_match.group(1).strip()
        sp_name = re.sub(r'[^a-zA-Z0-9_\s]', '', sp_name)
        return sp_name.replace(" ", "_")

    id_upper = record_id.upper()
    if any(id_upper.startswith(p) for p in ["PF3D7", "PF10", "PF11", "MAL", "PF"]):
        return "Plasmodium_falciparum"
    elif any(id_upper.startswith(p) for p in ["PVP01", "PVX", "PV"]):
        return "Plasmodium_vivax"
    elif id_upper.startswith("PBANKA") or id_upper.startswith("PB"):
        return "Plasmodium_berghei"
    elif id_upper.startswith("PCHAS") or id_upper.startswith("PCH"):
        return "Plasmodium_chabaudi"
    elif id_upper.startswith("PKNH") or id_upper.startswith("PK"):
        return "Plasmodium_knowlesi"
    elif id_upper.startswith("PY17X") or id_upper.startswith("PY"):
        return "Plasmodium_yoelii"
    elif id_upper.startswith("PCOAH"):
        return "Plasmodium_coatneyi"

    return "Plasmodium_sp"

def summarize_sequence(seq_str, gene_id, species, label):
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
    print("[ERROR] Could not locate 'promoter_sequences.fasta' or 'family_sequences.fasta'.")
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
            res = summarize_sequence(seq_str, rec.id.split()[0], species, label)
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

    for seq_type, group in master_table.groupby("Seq_Type"):
        df_melted = group.melt(
            id_vars=['Species', 'Gene_ID'],
            value_vars=property_cols,
            var_name='Property',
            value_name='Value'
        ).dropna()

        # Single Boxplot Canvas (hue='Property' added to suppress Seaborn warning)
        plt.figure(figsize=(14, 7), dpi=300)
        
        ax = sns.boxplot(
            data=df_melted, x='Property', y='Value', hue='Property',
            palette="Set3", showfliers=False, boxprops=dict(alpha=0.8),
            linewidth=1.2, legend=False
        )
        
        sns.stripplot(
            data=df_melted, x='Property', y='Value',
            color="black", alpha=0.3, size=3, jitter=0.25
        )

        plt.title(f"{seq_type} — Comparative DNA Structural Properties Profile", fontsize=14, fontweight='bold', pad=15)
        plt.xlabel("Structural Property", fontsize=12, fontweight='bold', labelpad=10)
        plt.ylabel("Value Range / Metric Score", fontsize=12, fontweight='bold', labelpad=10)
        plt.xticks(rotation=25, ha='right', fontsize=11)
        
        plt.tight_layout()
        plot_name = f"{seq_type}_ALL_PROPERTIES_UNIFIED_BOXPLOT.png"
        plt.savefig(os.path.join(PLOT_DIR, plot_name), dpi=300)
        plt.savefig(os.path.join(BASE_OUTPUT, plot_name), dpi=300)
        plt.close()

        # Faceted Grid Plot (hue='Species' added to suppress Seaborn warning)
        g = sns.catplot(
            data=df_melted, x='Species', y='Value', hue='Species', col='Property',
            col_wrap=3, kind='box', sharey=False, height=3.5, aspect=1.2,
            palette='Blues_d', showfliers=False, legend=False
        )
        g.set_xticklabels(rotation=45, ha='right', fontstyle='italic')
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle(f"{seq_type} — All DNA Structural Properties by Species", fontsize=14, fontweight='bold')
        grid_name = f"{seq_type}_ALL_PROPERTIES_FACETED_GRID.png"
        g.savefig(os.path.join(PLOT_DIR, grid_name), dpi=300)
        g.savefig(os.path.join(BASE_OUTPUT, grid_name), dpi=300)
        plt.close()

print("\n========================================")
print("Pipeline complete!")
print("Outputs yield folder created at: outputs/DNA_Structural_Property_Results/")
print("Plots generated in: outputs/DNA_Structural_Property_Results/Plots/ and outputs/")
print("========================================\n")