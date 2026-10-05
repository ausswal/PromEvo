#!/usr/bin/env python3
"""
NON-B DNA MOTIF ANALYSIS PIPELINE (WARNING-FREE PUBLICATION PLOTTING)
Authors: Swarup Das and Subarna Thakur
Affiliation: Department of Bioinformatics, University of North Bengal

Scans PROMOTER DNA only (family_sequences.fasta is protein and is NOT scanned:
A/C/G/T are amino-acid letters, and G{3,5} matches glycine runs, so any
"motif" counts from it were meaningless).

Counting units (both are "one hit = one non-overlapping feature"):
  TTS_count  = number of triplex target TRACTS (overlapping 15-bp windows that
               pass the test are merged into a single tract), both strands.
  QSFM_count = number of non-overlapping G-quadruplex matches, both strands.
  Total_NonB_Motifs = TTS_count + QSFM_count (a sum of like units).
  *_density  = the same counts per kb of sequence analysed (promoters are
               clipped at contig ends, so many are shorter than PROMOTER_LEN).
Promoters longer than PROMOTER_LEN are trimmed to their LAST PROMOTER_LEN bp
(TSS-anchored: the 3' end of each promoter FASTA entry sits at the TSS).
"""

import os
import glob
import re
import csv
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns

from pipeline_config import PROMOTER_LEN, check_promoter_lengths

# ==============================================================================
# 0. CONFIGURATION & DIRECTORY CREATION
# ==============================================================================
BASE_OUTPUT = "outputs"
OUT_ROOT = os.path.join(BASE_OUTPUT, "NonB_DNA_Motif_Results")
PLOT_DIR = os.path.join(OUT_ROOT, "Plots")

os.makedirs(BASE_OUTPUT, exist_ok=True)
os.makedirs(OUT_ROOT, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


def locate_sequence_files():
    """Locates the PROMOTER FASTA only (family_sequences.fasta is protein --
    see module docstring). Search order matches the rest of the pipeline
    (outputs/ first); the first hit wins so one file is never counted twice."""
    for r in ['outputs', '.', '..']:
        promoter_f = os.path.join(r, "promoter_sequences.fasta")
        if os.path.exists(promoter_f):
            return [("Promoters", promoter_f)]
    return []


# ==============================================================================
# 1. NON-B DNA DETECTION ALGORITHMS
# ==============================================================================

def detect_TTS(seq):
    """Triplex Target Site detection: 15-bp windows with >=50% G and <=1
    pyrimidine are found, then OVERLAPPING passing windows are merged into
    tracts. Returns the number of tracts, so a 30-bp qualifying tract counts
    as 1 (not 16 windows) -- the same 'one non-overlapping feature' unit as
    detect_QSFM."""
    seq = str(seq).upper()
    if len(seq) < 15:
        return 0
    tracts = 0
    cur_end = -1                      # end (exclusive) of the tract being built
    for i in range(len(seq) - 14):
        window = seq[i:i+15]
        g_count = window.count("G")
        pyrimidine_count = window.count("C") + window.count("T")
        if g_count / 15.0 >= 0.5 and pyrimidine_count <= 1:
            if i >= cur_end:          # not overlapping the current tract -> new tract
                tracts += 1
            cur_end = i + 15
    return tracts


def detect_QSFM(seq):
    """G-Quadruplex motif detection (Pattern: G3-5 N1-7 G3-5 N1-7 G3-5 N1-7 G3-5)"""
    seq = str(seq).upper()
    if len(seq) < 20:
        return 0
    pattern = re.compile(r"G{3,5}[ATCG]{1,7}G{3,5}[ATCG]{1,7}G{3,5}[ATCG]{1,7}G{3,5}")
    matches = pattern.findall(seq)
    return len(matches)


def scan_both_strands(seq, fn):
    """Scan both forward and reverse complement strands"""
    bio_seq = Seq(seq)
    hits_fwd = fn(str(bio_seq))
    hits_rev = fn(str(bio_seq.reverse_complement()))
    return hits_fwd + hits_rev


def load_organism_map(labels_path="outputs/family_labels.csv"):
    """Returns dict seq_id -> organism from family_labels.csv, the
    authoritative source already used throughout this pipeline. Preferred
    over extract_clean_species()'s FASTA-header heuristics below, which
    are Plasmodium-specific (its final fallback mislabels ANY unrecognized
    accession as 'Plasmodium_sp' -- silently wrong for Toxoplasma/Babesia
    sequences, whose accessions don't match any of the Plasmodium prefix
    checks). Returns {} if unavailable, so callers can fall back cleanly."""
    if not os.path.exists(labels_path):
        return {}
    try:
        df = pd.read_csv(labels_path)
    except Exception:
        return {}
    lower_cols = {c.lower(): c for c in df.columns}
    seq_id_col = lower_cols.get("seq_id")
    organism_col = next((lower_cols[c] for c in ("organism", "species", "organism_name") if c in lower_cols), None)
    if seq_id_col is None or organism_col is None:
        return {}
    return dict(zip(df[seq_id_col].astype(str), df[organism_col].astype(str)))


def extract_genus(species):
    """First token of the species string -- handles both underscore-joined
    ('Plasmodium_falciparum', this script's own convention) and
    space-joined ('Plasmodium falciparum', family_labels.csv's convention)
    forms."""
    species = (species or "").strip()
    if not species or species.lower() == "nan":
        return "Unknown_genus"
    return re.split(r"[_\s]+", species)[0]


def sanitize_filename(name):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip())
    return safe.strip("_") or "Unknown"


def extract_clean_species(record_id, record_desc):
    """
    Robust species parser specifically designed for Plasmodium FASTA headers.
    Prevents individual gene IDs from becoming x-axis ticks.
    """
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
    elif id_upper.startswith("PML"):
        return "Plasmodium_malariae"
    elif id_upper.startswith("PREL"):
        return "Plasmodium_relictum"

    plasmo_match = re.search(r'(Plasmodium_[a-z]+)', full_header, re.IGNORECASE)
    if plasmo_match:
        return plasmo_match.group(1)

    # No safe default exists: the old "Plasmodium_sp" mislabeled every
    # non-Plasmodium record whose ID matched none of the prefixes above.
    return "Unknown_species"


# ==============================================================================
# 2. MAIN PROCESSING LOOP
# ==============================================================================

targets = locate_sequence_files()

if not targets:
    print("[ERROR] Could not locate 'promoter_sequences.fasta' (looked in outputs/, ., ..).")
    exit(1)

organism_map = load_organism_map()
if organism_map:
    print(f"[i] Using authoritative organism mapping from family_labels.csv "
          f"({len(organism_map)} sequence(s) mapped).")
else:
    print("[!] No family_labels.csv organism mapping found — falling back to "
          "FASTA-header heuristics (Plasmodium-specific; non-Plasmodium accessions "
          "will fall back to 'Plasmodium_sp' — re-run with family_labels.csv present "
          "for correct Toxoplasma/Babesia labeling).")

n_heuristic_fallback = 0
all_results = []

for label, fpath in targets:
    print(f"\n==== Processing Non-B Motifs for Target: {label} ({os.path.basename(fpath)}) ====")
    
    try:
        records = list(SeqIO.parse(fpath, "fasta"))
        if not records:
            continue

        print(f"  • Found {len(records)} sequences")
        # trims=True: this stage trims longer promoters itself (below).
        check_promoter_lengths([len(r.seq) for r in records], stage="010", trims=True)

        for record in records:
            seq_str = str(record.seq).upper()
            if not seq_str:
                continue
            
            # TSS-anchored trim: keep the LAST PROMOTER_LEN bp (nearest the TSS).
            if len(seq_str) > PROMOTER_LEN:
                seq_str = seq_str[-PROMOTER_LEN:]

            gene_id = record.id.split()[0]
            species = organism_map.get(gene_id)
            if not species:
                species = extract_clean_species(record.id, record.description)
                n_heuristic_fallback += 1
            genus = extract_genus(species)

            tts_hits = scan_both_strands(seq_str, detect_TTS)
            qs_hits  = scan_both_strands(seq_str, detect_QSFM)
            
            seq_kb = len(seq_str) / 1000.0
            total_hits = tts_hits + qs_hits
            all_results.append({
                "Seq_Type": label,
                "Species": species,
                "Genus": genus,
                "Gene_ID": gene_id,
                "Sequence_Length": len(seq_str),
                "Sequence_Length_kb": seq_kb,
                "TTS_count": tts_hits,
                "QSFM_count": qs_hits,
                "Total_NonB_Motifs": total_hits,
                # hits per kb -- length-normalized, since clipped promoters
                # are shorter than PROMOTER_LEN
                "TTS_density": tts_hits / seq_kb,
                "QSFM_density": qs_hits / seq_kb,
                "Total_density": total_hits / seq_kb,
            })

    except Exception as e:
        print(f"  [!] ERROR processing {fpath}: {e}")

if n_heuristic_fallback:
    print(f"\n[!] {n_heuristic_fallback} sequence(s) had no family_labels.csv entry and "
          f"used the FASTA-header heuristic instead.")

if all_results:
    df_all = pd.DataFrame(all_results)
    
    master_csv = os.path.join(OUT_ROOT, "NonB_DNA_Motif_Counts_Master.csv")
    df_all.to_csv(master_csv, index=False)
    df_all.to_csv(os.path.join(BASE_OUTPUT, "nonB_dna_motif_counts.csv"), index=False)

    for seq_type, group in df_all.groupby("Seq_Type"):
        group.to_csv(os.path.join(OUT_ROOT, f"NonB_DNA_Motif_Counts_{seq_type}.csv"), index=False)

    # ==============================================================================
    # 3. HIGH-QUALITY PUBLICATION VISUALIZATIONS
    # ==============================================================================
    sns.set_theme(style="whitegrid", font="sans-serif")

    METRICS = [("TTS_count", "Triplex Target Sites (TTS)"),
               ("QSFM_count", "G-Quadruplex Motifs (QSFM)"),
               ("Total_NonB_Motifs", "Total Non-B DNA Motifs")]

    # ---- 3a. Cross-species comparison (unchanged from before -- still a
    # genuinely useful, different view: all species side-by-side in one plot) ----
    for seq_type, group in df_all.groupby("Seq_Type"):
        for metric_col, metric_label in METRICS:
            plt.figure(figsize=(10, 6), dpi=300)
            
            species_order = group.groupby('Species')[metric_col].mean().sort_values(ascending=False).index
            
            ax = sns.barplot(
                data=group, x='Species', y=metric_col, hue='Species', order=species_order,
                palette="Blues_d", errorbar='se', capsize=0.1, err_kws={'linewidth': 1.5}, legend=False
            )
            
            sns.stripplot(
                data=group, x='Species', y=metric_col, order=species_order,
                color="crimson", alpha=0.5, size=4, jitter=0.2
            )
            
            clean_species_labels = [s.replace("_", " ") for s in species_order]
            
            ax.set_xticks(range(len(clean_species_labels)))
            ax.set_xticklabels(clean_species_labels, rotation=30, ha='right', fontsize=10, fontstyle='italic')
            
            plt.title(f"{seq_type} — Mean {metric_label} per Gene Across Species", fontsize=12, fontweight='bold', pad=15)
            plt.xlabel("Species", fontsize=11, fontweight='bold', labelpad=10)
            plt.ylabel(f"Mean Count per Sequence (±SE)", fontsize=11, fontweight='bold', labelpad=10)
            
            plt.ylim(bottom=0)
            plt.tight_layout()
            
            plot_filename = f"{seq_type}_{metric_col}.png"
            plt.savefig(os.path.join(PLOT_DIR, plot_filename), dpi=300)
            plt.savefig(os.path.join(BASE_OUTPUT, f"nonB_{plot_filename}"), dpi=300)
            plt.close()

    # ---- 3b. Per-organism AND per-genus distribution plots (new) ----
    def plot_group_distribution(group_df, seq_type, label, metric_col, metric_label, out_dir, tag_word):
        """One figure showing the distribution of `metric_col` across every
        gene in this single group (one organism, or every species pooled
        under one genus) -- not a cross-group comparison, a within-group
        distribution, since that's what's actually informative once you've
        already split the data down to one organism/genus."""
        n_genes = len(group_df)
        plt.figure(figsize=(4, 6), dpi=200)
        ax = sns.boxplot(data=group_df, y=metric_col, color="#8FAADC", width=0.4,
                          showfliers=False, boxprops=dict(alpha=0.7))
        sns.stripplot(data=group_df, y=metric_col, color="crimson", alpha=0.6, size=4, jitter=0.15, ax=ax)
        mean_val = group_df[metric_col].mean()
        ax.axhline(mean_val, color="black", linestyle="--", linewidth=1, label=f"mean={mean_val:.2f}")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylim(bottom=0)
        ax.set_xlabel("")
        ax.set_xticks([])
        ax.set_ylabel(f"{metric_label} per gene", fontsize=10, fontweight="bold")
        clean_label = label.replace("_", " ")
        plt.title(f"{seq_type} — {tag_word.capitalize()}: {clean_label}\n(n={n_genes} gene(s))",
                  fontsize=10.5, fontweight="bold")
        plt.tight_layout()
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{sanitize_filename(label)}_{metric_col}.png")
        plt.savefig(out_path, dpi=200)
        plt.close()
        return out_path

    print(f"\n{'#' * 60}\n# PER-ORGANISM DISTRIBUTION PLOTS\n{'#' * 60}")
    for seq_type, seq_group in df_all.groupby("Seq_Type"):
        org_out_dir = os.path.join(OUT_ROOT, "By_Organism", seq_type)
        for species, sp_group in seq_group.groupby("Species"):
            for metric_col, metric_label in METRICS:
                out_path = plot_group_distribution(sp_group, seq_type, species, metric_col,
                                                    metric_label, org_out_dir, tag_word="organism")
            print(f"  [✓] {seq_type} / {species}: {len(sp_group)} gene(s) -> {org_out_dir}/")

    print(f"\n{'#' * 60}\n# PER-GENUS DISTRIBUTION PLOTS (pooled across species)\n{'#' * 60}")
    for seq_type, seq_group in df_all.groupby("Seq_Type"):
        genus_out_dir = os.path.join(OUT_ROOT, "By_Genus", seq_type)
        for genus, genus_group in seq_group.groupby("Genus"):
            n_species_in_genus = genus_group["Species"].nunique()
            for metric_col, metric_label in METRICS:
                out_path = plot_group_distribution(genus_group, seq_type, genus, metric_col,
                                                    metric_label, genus_out_dir, tag_word="genus")
            print(f"  [✓] {seq_type} / {genus}: {len(genus_group)} gene(s) across "
                  f"{n_species_in_genus} species -> {genus_out_dir}/")

    # ==============================================================================
    # 4. ADDITIONAL VISUALIZATIONS
    # ==============================================================================
    # All six use only columns already present in df_all -- no change to the
    # detection logic (detect_TTS / detect_QSFM) or to any plot above this point.

    # (Sequence_Length_kb and the *_density columns are computed per sequence
    # above, so every CSV written earlier already contains them.)

    DENSITY_METRICS = [("TTS_density", "TTS per kb"),
                        ("QSFM_density", "QSFM per kb"),
                        ("Total_density", "Total Non-B Motifs per kb")]

    def plot_cross_species_generic(group, seq_type, metric_col, metric_label, y_label,
                                    out_dir, file_suffix, palette="Greens_d"):
        """Bar+strip cross-species plot, structurally identical to section 3a but
        parametrized so it can be reused for density metrics without touching 3a."""
        plt.figure(figsize=(10, 6), dpi=300)
        order = group.groupby('Species')[metric_col].mean().sort_values(ascending=False).index
        ax = sns.barplot(data=group, x='Species', y=metric_col, hue='Species', order=order,
                          palette=palette, errorbar='se', capsize=0.1,
                          err_kws={'linewidth': 1.5}, legend=False)
        sns.stripplot(data=group, x='Species', y=metric_col, order=order,
                       color="crimson", alpha=0.5, size=4, jitter=0.2)
        clean_labels = [s.replace("_", " ") for s in order]
        ax.set_xticks(range(len(clean_labels)))
        ax.set_xticklabels(clean_labels, rotation=30, ha='right', fontsize=10, fontstyle='italic')
        plt.title(f"{seq_type} — Mean {metric_label} per Gene Across Species",
                  fontsize=12, fontweight='bold', pad=15)
        plt.xlabel("Species", fontsize=11, fontweight='bold', labelpad=10)
        plt.ylabel(y_label, fontsize=11, fontweight='bold', labelpad=10)
        plt.ylim(bottom=0)
        plt.tight_layout()
        os.makedirs(out_dir, exist_ok=True)
        fname = f"{seq_type}_{file_suffix}.png"
        plt.savefig(os.path.join(out_dir, fname), dpi=300)
        plt.savefig(os.path.join(BASE_OUTPUT, f"nonB_{fname}"), dpi=300)
        plt.close()

    # ---- 4a. Motif density per kb (length-normalized) ----
    density_dir = os.path.join(OUT_ROOT, "Density")
    for seq_type, group in df_all.groupby("Seq_Type"):
        for metric_col, metric_label in DENSITY_METRICS:
            plot_cross_species_generic(group, seq_type, metric_col, metric_label,
                                        f"Mean {metric_label} (±SE)", density_dir, metric_col)
    print(f"[✓] Density plots (per kb, length-normalized) -> {density_dir}/")

    # ---- 4b. TTS vs QSFM composition (stacked bar) per species ----
    comp_dir = os.path.join(OUT_ROOT, "Composition")
    os.makedirs(comp_dir, exist_ok=True)
    for seq_type, group in df_all.groupby("Seq_Type"):
        means = group.groupby('Species')[['TTS_count', 'QSFM_count']].mean()
        means = means.loc[means.sum(axis=1).sort_values(ascending=False).index]
        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
        x = np.arange(len(means))
        ax.bar(x, means['TTS_count'], label='TTS', color='#4C72B0')
        ax.bar(x, means['QSFM_count'], bottom=means['TTS_count'], label='QSFM', color='#DD8452')
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace('_', ' ') for s in means.index], rotation=30,
                            ha='right', fontstyle='italic', fontsize=10)
        ax.set_ylabel("Mean Count per Sequence", fontsize=11, fontweight='bold')
        ax.set_xlabel("Species", fontsize=11, fontweight='bold')
        ax.set_title(f"{seq_type} — TTS vs QSFM Composition by Species",
                     fontsize=12, fontweight='bold', pad=15)
        ax.legend(title="Motif Type", fontsize=9)
        plt.tight_layout()
        fname = f"{seq_type}_TTS_QSFM_composition.png"
        plt.savefig(os.path.join(comp_dir, fname), dpi=300)
        plt.savefig(os.path.join(BASE_OUTPUT, f"nonB_{fname}"), dpi=300)
        plt.close()
    print(f"[✓] TTS vs QSFM composition plots -> {comp_dir}/")

    # ---- 4d. Sequence length vs motif count scatter (with Pearson r) ----
    scatter_dir = os.path.join(OUT_ROOT, "Length_Scatter")
    os.makedirs(scatter_dir, exist_ok=True)
    for seq_type, group in df_all.groupby("Seq_Type"):
        for metric_col, metric_label in METRICS:
            plt.figure(figsize=(7.5, 6), dpi=300)
            ax = sns.scatterplot(data=group, x='Sequence_Length', y=metric_col, hue='Genus',
                                  palette='tab10', s=55, alpha=0.85, edgecolor='black', linewidth=0.3)
            if group['Sequence_Length'].nunique() > 1 and group[metric_col].std() > 0:
                sns.regplot(data=group, x='Sequence_Length', y=metric_col, scatter=False,
                            color='black', line_kws={'linewidth': 1.5, 'linestyle': '--'}, ax=ax)
                r = np.corrcoef(group['Sequence_Length'], group[metric_col])[0, 1]
                ax.text(0.02, 0.97, f"Pearson r = {r:.2f}", transform=ax.transAxes,
                        fontsize=9, va='top', bbox=dict(boxstyle='round', fc='white', alpha=0.85))
            plt.title(f"{seq_type} — {metric_label} vs Sequence Length",
                      fontsize=12, fontweight='bold', pad=12)
            plt.xlabel("Sequence Length (bp)", fontsize=11, fontweight='bold')
            plt.ylabel(metric_label, fontsize=11, fontweight='bold')
            plt.legend(title="Genus", fontsize=8, bbox_to_anchor=(1.02, 1), loc='upper left')
            plt.tight_layout()
            fname = f"{seq_type}_{metric_col}_vs_length.png"
            plt.savefig(os.path.join(scatter_dir, fname), dpi=300)
            plt.savefig(os.path.join(BASE_OUTPUT, f"nonB_{fname}"), dpi=300)
            plt.close()
    print(f"[✓] Sequence-length-vs-motif-count scatter plots -> {scatter_dir}/")

    # ---- 4e. Species/genus x metric heatmap ----
    heatmap_dir = os.path.join(OUT_ROOT, "Heatmaps")
    os.makedirs(heatmap_dir, exist_ok=True)
    for seq_type, group in df_all.groupby("Seq_Type"):
        for level_col in ("Species", "Genus"):
            pivot = group.groupby(level_col)[['TTS_count', 'QSFM_count', 'Total_NonB_Motifs']].mean()
            pivot = pivot.sort_values('Total_NonB_Motifs', ascending=False)
            pivot.index = [i.replace('_', ' ') for i in pivot.index]
            plt.figure(figsize=(6.5, max(3, 0.55 * len(pivot) + 1.2)), dpi=300)
            sns.heatmap(pivot, annot=True, fmt=".2f", cmap="YlGnBu", linewidths=0.5,
                        linecolor='white', cbar_kws={'label': 'Mean Count per Sequence'})
            plt.title(f"{seq_type} — Mean Non-B Motif Counts by {level_col}",
                      fontsize=12, fontweight='bold', pad=12)
            plt.ylabel(level_col, fontsize=11, fontweight='bold')
            plt.xlabel("Metric", fontsize=11, fontweight='bold')
            plt.xticks(rotation=20, ha='right')
            plt.yticks(rotation=0, fontstyle='italic' if level_col == 'Species' else 'normal')
            plt.tight_layout()
            fname = f"{seq_type}_{level_col}_heatmap.png"
            plt.savefig(os.path.join(heatmap_dir, fname), dpi=300)
            plt.savefig(os.path.join(BASE_OUTPUT, f"nonB_{fname}"), dpi=300)
            plt.close()
    print(f"[✓] Species/genus x metric heatmaps -> {heatmap_dir}/")

    # ---- 4f. Violin plots (distribution shape) cross-species ----
    violin_dir = os.path.join(OUT_ROOT, "Violin")
    os.makedirs(violin_dir, exist_ok=True)
    for seq_type, group in df_all.groupby("Seq_Type"):
        for metric_col, metric_label in METRICS:
            plt.figure(figsize=(10, 6), dpi=300)
            order = group.groupby('Species')[metric_col].mean().sort_values(ascending=False).index
            ax = sns.violinplot(data=group, x='Species', y=metric_col, hue='Species', order=order,
                                 palette="Purples", inner="quartile", cut=0, legend=False)
            sns.stripplot(data=group, x='Species', y=metric_col, order=order,
                          color="black", alpha=0.4, size=3, jitter=0.15)
            clean_labels = [s.replace('_', ' ') for s in order]
            ax.set_xticks(range(len(clean_labels)))
            ax.set_xticklabels(clean_labels, rotation=30, ha='right', fontstyle='italic', fontsize=10)
            plt.title(f"{seq_type} — {metric_label} Distribution by Species",
                      fontsize=12, fontweight='bold', pad=15)
            plt.xlabel("Species", fontsize=11, fontweight='bold')
            plt.ylabel(metric_label, fontsize=11, fontweight='bold')
            plt.ylim(bottom=0)
            plt.tight_layout()
            fname = f"{seq_type}_{metric_col}_violin.png"
            plt.savefig(os.path.join(violin_dir, fname), dpi=300)
            plt.savefig(os.path.join(BASE_OUTPUT, f"nonB_{fname}"), dpi=300)
            plt.close()
    print(f"[✓] Violin distribution plots -> {violin_dir}/")

print("\n========================================")
print("Pipeline complete!")
print("Outputs folder created at: outputs/NonB_DNA_Motif_Results/")
print("Plots saved in: outputs/NonB_DNA_Motif_Results/Plots/ and outputs/")
print("========================================\n")