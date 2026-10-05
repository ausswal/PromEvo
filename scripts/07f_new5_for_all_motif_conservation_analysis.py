#!/usr/bin/env python3
"""
STAGE 14 (v4) — Cross-Species / Cross-Paralog Motif Conservation Analysis,
now run per-ORGANISM and per-GENUS (matching 07c's per-organism/per-genus
meme.xml split and the same discovery pattern used in
07d_new_alt_fimo_goatools.py).

WHAT CHANGED FROM THE PREVIOUS (single combined) VERSION
-----------------------------------------------------------
07c now writes one meme.xml per organism (outputs/meme_output/<Organism>/)
and one per genus (outputs/meme_output_genus/<Genus>/), instead of a
single outputs/meme_output/meme.xml. This auto-discovers every meme.xml
under BOTH locations (glob, no hardcoded organism/genus names -- same
discover_meme_xml_groups() approach as 07d) and runs the full existing
conservation-analysis pipeline (parsing, species mapping, MCS scoring,
Kruskal-Wallis, the four-panel figure) independently for each, producing:

  outputs/motif_conservation_by_organism/<Organism>_motif_conservation_analysis.png
  outputs/motif_conservation_by_genus/<Genus>_motif_conservation_analysis.png
  outputs/motif_conservation_summary_by_organism.csv   (combined, tagged `organism`)
  outputs/motif_conservation_summary_by_genus.csv      (combined, tagged `genus`)
  outputs/motif_occurrences_by_organism.csv            (combined, tagged `organism`)
  outputs/motif_occurrences_by_genus.csv               (combined, tagged `genus`)

07c_new1_Motif_discovery_for_all.py writes exactly the two directory
layouts above (verified against the 07c source). If neither exists, this
falls back to the single legacy outputs/meme_output/meme.xml.

COORDINATE SYSTEM (v4 FIX) -- TSS-RELATIVE POSITIONS
------------------------------------------------------
MEME reports each site's `position` as a 1-based offset from the START of
the (upstream-only) promoter record, i.e. position 1 = the base furthest
upstream (-L) and position L = the base immediately before the TSS / gene
start (-1). Using that raw number for "positional constraint" is wrong in
two ways: it reads backwards biologically (a larger number is CLOSER to the
TSS), and it mis-aligns promoters of different length (promoters clipped at
a contig edge are shorter than 1000 bp, so the same raw number means a
different distance from the TSS).

This version converts every site to TSS-relative coordinates using each
sequence's own length from meme.xml's <training_set>:

    pos_rel_tss_start  = position - seq_len - 1      (last base -> -1)
    pos_rel_tss_end    = pos_rel_tss_start + width - 1
    pos_rel_tss_center = pos_rel_tss_start + (width - 1) / 2

All positional statistics (positional_std_bp, MCS, Kruskal-Wallis) and the
Panel B plot use pos_rel_tss_start. The raw MEME value is kept in the
output CSV as `position` (unchanged, for backward compatibility) next to
the new `pos_rel_tss_*` columns.

Assumption to spot-check once: for sites MEME reports on the minus strand
(only present if 07c runs MEME with -revcomp), `position` is taken to be
the leftmost base of the site on the forward strand as given in the FASTA.

IMPORTANT REFRAMING FOR ORGANISM-LEVEL RESULTS
------------------------------------------------
The original analysis' whole point was measuring how conserved a
MEME-discovered motif is ACROSS SPECIES. That framing only holds at the
GENUS level now, where a genus's meme.xml was built from multiple
species' sequences. An ORGANISM-level meme.xml was built from just ONE
organism's own sequences (its paralogous gene-family members, if it has
more than one copy) -- so an organism-level plot is really showing
CROSS-PARALOG conservation WITHIN one genome, not cross-species
conservation. The "species" axis on Panels A/D will just show one
species name repeated across that organism's own paralogs. This isn't a
bug -- it's a real difference in what's being measured, and it's called
out in each organism-level plot's title/panel labels and in the console
output so it's never confused with the genus-level (genuinely
cross-species) result.

Organisms with too few sequences for MEME to have discovered a motif at
all (most commonly single-copy genes) are skipped with a clear message,
not silently dropped.

BACKWARD-COMPATIBILITY NOTE: the previous single-run version wrote to
outputs/motif_conservation_summary.csv, outputs/motif_occurrences_per_species.csv,
and outputs/motif_conservation_analysis.png. This version does NOT write
those exact paths -- if 07h_motif_turnover_analysis.py or anything else
reads those specific filenames, it will need pointing at the new
per-organism/per-genus outputs instead.

Everything else (MEME XML parsing, species-mapping fix, simulated-data
guard, all four plot panels, layout fixes) is unchanged from the previous
corrected version -- just re-run once per discovered group instead of
once globally.
"""

import os
import re
import sys
import glob
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from scipy import stats

# ============================== CONFIGURATION ==============================
MEME_OUT_DIR = "outputs/meme_output"                  # per-organism (07c)
MEME_OUT_DIR_GENUS = "outputs/meme_output_genus"       # per-genus (07c)
LEGACY_MEME_XML = "outputs/meme_output/meme.xml"       # pre-split fallback

FAMILY_LABELS_CSV = "outputs/family_labels.csv"

PLOT_DIR = "outputs/motif_conservation_by_organism"
PLOT_DIR_GENUS = "outputs/motif_conservation_by_genus"

OUT_SUMMARY_CSV = "outputs/motif_conservation_summary_by_organism.csv"
OUT_SUMMARY_CSV_GENUS = "outputs/motif_conservation_summary_by_genus.csv"
OUT_OCCURRENCE_CSV = "outputs/motif_occurrences_by_organism.csv"
OUT_OCCURRENCE_CSV_GENUS = "outputs/motif_occurrences_by_genus.csv"

ALLOW_SIMULATED_FALLBACK = False

# Used ONLY if a meme.xml <sequence> entry has no usable length attribute.
PROMOTER_LENGTH_DEFAULT = 1000
# ===========================================================================


def sanitize_filename(name):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip())
    return safe.strip("_") or "Unknown"


def discover_meme_xml_groups(meme_out_dir):
    """Finds every <meme_out_dir>/<Label>/meme.xml and returns {label: path}.
    Same discovery approach as 07d_new_alt_fimo_goatools.py -- no
    hardcoded organism/genus names."""
    groups = {}
    for xml_path in sorted(glob.glob(os.path.join(meme_out_dir, "*", "meme.xml"))):
        label = os.path.basename(os.path.dirname(xml_path))
        groups[label] = xml_path
    return groups


def parse_meme_xml(xml_path):
    """Unchanged from the previous version."""
    if not os.path.exists(xml_path):
        return None, None
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()

        seq_map = {}
        seq_len_map = {}
        for seq in root.findall(".//training_set/sequence"):
            seq_map[seq.get("id")] = seq.get("name")
            try:
                seq_len_map[seq.get("id")] = int(seq.get("length"))
            except (TypeError, ValueError):
                seq_len_map[seq.get("id")] = None

        motifs_data, occurrences = [], []
        for motif in root.findall(".//motifs/motif"):
            m_id = motif.get("id")
            m_name = motif.get("name")
            m_width = int(motif.get("width"))
            motifs_data.append({
                "motif_id": m_id,
                "motif_name": m_name,
                "width": m_width,
                "sites_count": int(motif.get("sites")),
                "e_value": float(motif.get("e_value")),
                "consensus": motif.get("alt") or m_name
            })
            for site in motif.findall(".//contributing_sites/contributing_site"):
                seq_id = site.get("sequence_id")
                occurrences.append({
                    "motif_id": m_id,
                    "motif_name": m_name,
                    "seq_id": seq_map.get(seq_id, seq_id),
                    "position": int(site.get("position")),   # raw MEME, 1-based from record start
                    "width": m_width,
                    "seq_len": seq_len_map.get(seq_id),
                    "pvalue": float(site.get("pvalue")),
                    "strand": site.get("strand", "+")
                })

        df_motifs = pd.DataFrame(motifs_data)
        df_occ = pd.DataFrame(occurrences)
        df_motifs["data_source"] = "REAL_MEME"
        df_occ["data_source"] = "REAL_MEME"
        return df_motifs, df_occ
    except Exception as e:
        print(f"    Warning: Failed to parse {xml_path}: {e}")
        return None, None


def generate_simulated_meme_data():
    """Unchanged."""
    print("=" * 70)
    print("[!!!] SIMULATED DEMO DATA IN USE -- THESE ARE NOT REAL RESULTS [!!!]")
    print("=" * 70)
    motifs_data = [
        {"motif_id": "motif_1", "motif_name": "Motif 1 (AP2-Domain)", "width": 12, "sites_count": 28, "e_value": 1.2e-15, "consensus": "TGCATGCA"},
        {"motif_id": "motif_2", "motif_name": "Motif 2 (Promoter-Box)", "width": 10, "sites_count": 22, "e_value": 3.4e-08, "consensus": "TATAAAAG"},
        {"motif_id": "motif_3", "motif_name": "Motif 3 (Upstream Element)", "width": 8, "sites_count": 15, "e_value": 5.1e-04, "consensus": "CCGGCCGG"},
    ]
    species_list = [f"Species_{i}" for i in range(1, 10)]
    occurrences = []
    np.random.seed(42)
    for m in motifs_data:
        for sp in species_list:
            if np.random.rand() > 0.2:
                occurrences.append({
                    "motif_id": m["motif_id"], "motif_name": m["motif_name"],
                    "seq_id": f"{sp}_gene_1", "species": sp,
                    "position": int(np.random.normal(loc=150, scale=30)),
                    "width": m["width"],
                    "seq_len": PROMOTER_LENGTH_DEFAULT,
                    "pvalue": float(np.random.uniform(1e-6, 1e-3)),
                    "strand": "+" if np.random.rand() > 0.5 else "-"
                })
    df_motifs = pd.DataFrame(motifs_data)
    df_occ = pd.DataFrame(occurrences)
    df_motifs["data_source"] = "SIMULATED_DEMO"
    df_occ["data_source"] = "SIMULATED_DEMO"
    return df_motifs, df_occ


def add_tss_relative_coordinates(df_occ, default_len=PROMOTER_LENGTH_DEFAULT):
    """Converts MEME's raw 1-based record-start positions to TSS-relative
    coordinates (last base of the upstream-only promoter = -1).

        pos_rel_tss_start  = position - seq_len - 1
        pos_rel_tss_end    = pos_rel_tss_start + width - 1
        pos_rel_tss_center = pos_rel_tss_start + (width - 1) / 2

    e.g. a 10-bp site starting at position 991 of a 1000-bp record covers
    the last 10 bases: start -10, end -1.
    """
    df = df_occ.copy()
    if "seq_len" not in df.columns:
        df["seq_len"] = np.nan
    missing = df["seq_len"].isna()
    if missing.any():
        print(f"    [!] {int(missing.sum())} site(s) had no sequence length in meme.xml -- "
              f"assuming {default_len} bp for those.")
        df.loc[missing, "seq_len"] = default_len
    df["seq_len"] = df["seq_len"].astype(int)

    df["pos_rel_tss_start"] = df["position"] - df["seq_len"] - 1
    df["pos_rel_tss_end"] = df["pos_rel_tss_start"] + df["width"] - 1
    df["pos_rel_tss_center"] = df["pos_rel_tss_start"] + (df["width"] - 1) / 2.0

    # Sanity check: every site must lie fully inside [-seq_len, -1].
    bad = (df["pos_rel_tss_end"] > -1) | (df["pos_rel_tss_start"] < -df["seq_len"])
    if bad.any():
        print(f"    [!] {int(bad.sum())} site(s) fall outside the expected [-seq_len, -1] "
              f"range after conversion -- check that meme.xml was built from the "
              f"upstream-only promoter FASTA.")
    return df


def build_species_map(df_occ, labels_path):
    """Unchanged."""
    if "species" in df_occ.columns:
        return df_occ
    if not os.path.exists(labels_path):
        print(f"    [!] {labels_path} not found. Treating each seq_id as its own species bucket.")
        df_occ = df_occ.copy()
        df_occ["species"] = df_occ["seq_id"]
        return df_occ

    labels_df = pd.read_csv(labels_path)
    lower_cols = {c.lower(): c for c in labels_df.columns}
    seq_id_col = lower_cols.get("seq_id")
    organism_col = next((lower_cols[c] for c in ("organism", "species", "organism_name") if c in lower_cols), None)

    if seq_id_col is None or organism_col is None:
        print(f"    [!] {labels_path} missing usable seq_id/organism columns. "
              f"Falling back to seq_id as species bucket.")
        df_occ = df_occ.copy()
        df_occ["species"] = df_occ["seq_id"]
        return df_occ

    spec_map = dict(zip(labels_df[seq_id_col], labels_df[organism_col]))
    df_occ = df_occ.copy()
    df_occ["species"] = df_occ["seq_id"].map(spec_map)
    mapped = df_occ["species"].notna().sum()
    total = len(df_occ)
    print(f"    [*] Species mapping: {mapped}/{total} sequences mapped from {labels_path}.")
    if mapped < total:
        mask = df_occ["species"].isna()
        df_occ.loc[mask, "species"] = "Unmapped_" + df_occ.loc[mask, "seq_id"].astype(str)
    return df_occ


def run_conservation_analysis_for_group(label, meme_xml_path, plot_out_dir, tag_word):
    """
    Runs the FULL conservation-analysis pipeline (parse, species-map, MCS
    scoring, Kruskal-Wallis, four-panel figure) for ONE group (one
    organism's or one genus's meme.xml). Returns (summary_df, occ_df),
    both tagged with a `tag_word` column set to `label`, or (None, None)
    if this group couldn't be analyzed.
    """
    print(f"\n{'-' * 60}\n{tag_word.capitalize()}: {label}\n{'-' * 60}")

    df_motifs, df_occ = parse_meme_xml(meme_xml_path)
    if df_motifs is None or df_motifs.empty:
        if not ALLOW_SIMULATED_FALLBACK:
            print(f"    [!] No real MEME motifs found at '{meme_xml_path}' -- skipping "
                  f"{tag_word} '{label}' (most likely too few sequences for MEME to have "
                  f"found a motif here).")
            return None, None
        df_motifs, df_occ = generate_simulated_meme_data()

    # SEPARATE from the df_motifs-empty case above: MEME can report at
    # least one motif but zero contributing_sites for it (seen in practice
    # on the Toxoplasma genus run -- a real MEME outcome for a group with
    # very few/short sequences, not a parsing bug). pd.DataFrame([]) then
    # has NO COLUMNS AT ALL, so df_occ["seq_id"] below would raise
    # KeyError('seq_id') before ever reaching .map() -- this is exactly
    # what crashed the pipeline. Treat it the same way: skip cleanly with
    # a clear message rather than crash.
    if df_occ is None or df_occ.empty or "seq_id" not in df_occ.columns:
        if not ALLOW_SIMULATED_FALLBACK:
            print(f"    [!] MEME found {len(df_motifs)} motif(s) in '{meme_xml_path}' but "
                  f"reported zero contributing sites for any of them -- skipping {tag_word} "
                  f"'{label}' (nothing to compute conservation from).")
            return None, None
        df_motifs, df_occ = generate_simulated_meme_data()

    using_simulated = (df_motifs["data_source"] == "SIMULATED_DEMO").any()
    df_occ = add_tss_relative_coordinates(df_occ)
    df_occ = build_species_map(df_occ, FAMILY_LABELS_CSV)
    total_species_count = df_occ["species"].nunique()
    is_single_species_organism = (tag_word == "organism" and total_species_count <= 1)

    if is_single_species_organism:
        print(f"    [i] Only 1 species present in this organism-level meme.xml -- this plot "
              f"shows CROSS-PARALOG conservation within '{label}''s own genome, not "
              f"cross-species conservation (see module docstring).")

    summary_list = []
    for m_id, group in df_occ.groupby("motif_id"):
        m_info = df_motifs[df_motifs["motif_id"] == m_id].iloc[0]
        species_present = group["species"].nunique()
        species_penetrance = species_present / max(1, total_species_count)
        pos_std = group["pos_rel_tss_start"].std() if len(group) > 1 else 0.0
        pos_std = 0.0 if pd.isna(pos_std) else pos_std
        mcs = (species_penetrance * (1.0 / (1.0 + (pos_std / 50.0)))) * 100.0
        summary_list.append({
            "motif_id": m_id, "motif_name": m_info["motif_name"],
            "consensus": m_info["consensus"], "e_value": m_info["e_value"],
            "total_instances": len(group), "species_present_count": species_present,
            "total_species_analyzed": total_species_count,
            "species_penetrance_pct": round(species_penetrance * 100, 2),
            "positional_std_bp": round(pos_std, 2),
            "mean_pos_rel_tss_bp": round(group["pos_rel_tss_start"].mean(), 1),
            "median_pos_rel_tss_bp": round(group["pos_rel_tss_start"].median(), 1),
            "motif_conservation_score": round(mcs, 2),
            "conservation_status": "Highly Conserved" if mcs >= 65 else ("Moderately Conserved" if mcs >= 35 else "Variable / Lineage Specific"),
            "data_source": "SIMULATED_DEMO" if using_simulated else "REAL_MEME"
        })
    summary_df = pd.DataFrame(summary_list).sort_values("motif_conservation_score", ascending=False)

    motif_pos_groups = [g["pos_rel_tss_start"].values for _, g in df_occ.groupby("motif_id") if len(g) > 1]
    kw_stat, p_val = stats.kruskal(*motif_pos_groups) if len(motif_pos_groups) > 1 else (0.0, 1.0)

    # ==================== four-panel figure (layout unchanged) ====================
    sns.set_theme(style="whitegrid", font="sans-serif")
    fig = plt.figure(figsize=(24, 7), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.75, 0.85, 1.35])

    species_order = sorted(df_occ["species"].unique())
    palette = sns.color_palette("Set2", n_colors=df_occ["motif_id"].nunique())
    motif_color_map = dict(zip(sorted(df_occ["motif_id"].unique()), palette))

    # --- PANEL B ---
    ax2 = fig.add_subplot(gs[0])
    sns.barplot(data=summary_df, x="motif_id", y="motif_conservation_score",
                hue="motif_id", palette=motif_color_map, legend=False,
                edgecolor="black", linewidth=0.8, ax=ax2)
    ax2.set_title("A. Motif Conservation Score (MCS)", fontsize=11, fontweight="bold", pad=10)
    ax2.set_xlabel("Motif ID", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Conservation Score (0-100)", fontsize=10, fontweight="bold")
    ax2.set_ylim(0, 115)
    ax2.tick_params(axis='x', rotation=30, labelsize=8)
    for lbl in ax2.get_xticklabels():
        lbl.set_ha('right')
    for p in ax2.patches:
        val = p.get_height()
        if val > 0:
            ax2.annotate(f"{val:.1f}", (p.get_x() + p.get_width() / 2., val + 2),
                         ha='center', va='bottom', fontsize=8, fontweight='bold')

    # --- PANEL C ---
    ax3 = fig.add_subplot(gs[1])
    sns.boxplot(data=df_occ, x="motif_id", y="pos_rel_tss_start", hue="motif_id",
                palette=motif_color_map, width=0.4, ax=ax3, legend=False,
                boxprops=dict(alpha=0.8))
    sns.stripplot(data=df_occ, x="motif_id", y="pos_rel_tss_start", color="black",
                  alpha=0.3, jitter=0.15, size=3.5, ax=ax3)
    ax3.set_title("B. Positional Constraint & Jitter", fontsize=11, fontweight="bold", pad=10)
    ax3.set_xlabel("Motif ID", fontsize=10, fontweight="bold")
    ax3.set_ylabel("Site start relative to TSS (bp; -1 = base just upstream)",
                   fontsize=10, fontweight="bold")
    max_len = int(df_occ["seq_len"].max())
    ax3.set_ylim(-max_len - 20, 20)
    ax3.axhline(0, color="red", linestyle=":", linewidth=1.3)
    ax3.text(0.99, 0.985, "TSS (0)", transform=ax3.transAxes, ha="right", va="top",
             fontsize=7.5, color="red")
    ax3.tick_params(axis='x', rotation=30, labelsize=8)
    for lbl in ax3.get_xticklabels():
        lbl.set_ha('right')
    p_text = f"p = {p_val:.2e}" if p_val < 0.001 else f"p = {p_val:.4f}"
    ax3.text(0.5, 1.06, f"Kruskal-Wallis: {p_text}", transform=ax3.transAxes,
              ha="center", va="bottom", fontsize=7.5,
              bbox=dict(boxstyle="round,pad=0.25", facecolor="#E5E9F0", edgecolor="black", alpha=0.85))

    # --- PANEL D ---
    ax4 = fig.add_subplot(gs[2])
    pivot = (
        df_occ.groupby(["species", "motif_id"]).size().unstack(fill_value=0)
        .reindex(index=species_order)
        .reindex(columns=sorted(df_occ["motif_id"].unique()), fill_value=0)
    )
    sns.heatmap(pivot, ax=ax4, cmap="YlGnBu", annot=True, fmt="d",
                annot_kws={"size": 7}, cbar_kws={"label": "Motif instances", "shrink": 0.8},
                linewidths=0.5, linecolor="white")
    panel_d_title = ("C. Paralog x Motif Heatmap" if is_single_species_organism
                      else "C. Species x Motif Conservation Heatmap")
    ax4.set_title(panel_d_title, fontsize=11, fontweight="bold", pad=10)
    ax4.set_xlabel("Motif ID", fontsize=10, fontweight="bold")
    ax4.set_ylabel("Species", fontsize=10, fontweight="bold")
    ax4.tick_params(axis='x', labelsize=8, rotation=30)
    ax4.tick_params(axis='y', labelsize=7, rotation=0)
    for lbl in ax4.get_xticklabels():
        lbl.set_ha('right')

    title_suffix = "  [!!! SIMULATED DEMO DATA -- NOT REAL RESULTS !!!]" if using_simulated else ""
    reframe_suffix = "  [cross-PARALOG, single species]" if is_single_species_organism else ""
    fig.suptitle(
        f"{tag_word.capitalize()}: {label} \u2014 Motif Conservation & Positional Constraint"
        f"{reframe_suffix}{title_suffix}",
        fontsize=13, fontweight="bold", color="red" if using_simulated else "black"
    )

    os.makedirs(plot_out_dir, exist_ok=True)
    slug = sanitize_filename(label)
    plot_path = os.path.join(plot_out_dir, f"{slug}_motif_conservation_analysis.png")
    plt.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"    [\u2713] Figure saved -> {plot_path}")

    for _, r in summary_df.iterrows():
        print(f"      - {r['motif_id']} ({r['consensus']}): MCS = {r['motif_conservation_score']} "
              f"| Status: {r['conservation_status']}")

    summary_df = summary_df.copy()
    summary_df[tag_word] = label
    df_occ_tagged = df_occ.copy()
    df_occ_tagged[tag_word] = label
    return summary_df, df_occ_tagged


def main():
    os.makedirs("outputs", exist_ok=True)

    species_groups = discover_meme_xml_groups(MEME_OUT_DIR)
    genus_groups = discover_meme_xml_groups(MEME_OUT_DIR_GENUS)

    if not species_groups and not genus_groups:
        if os.path.exists(LEGACY_MEME_XML):
            print(f"[!] No per-organism/per-genus meme.xml files found under {MEME_OUT_DIR}*/ -- "
                  f"falling back to the single legacy file at {LEGACY_MEME_XML}.")
            species_groups = {"All_organisms": LEGACY_MEME_XML}
        else:
            sys.exit(f"[!] No meme.xml files found under '{MEME_OUT_DIR}/<Organism>/', "
                      f"'{MEME_OUT_DIR_GENUS}/<Genus>/', or '{LEGACY_MEME_XML}'. "
                      f"Run 07c_Motif_discovery_for_all.py first.")

    print(f"[+] Found {len(species_groups)} organism-level and {len(genus_groups)} "
          f"genus-level meme.xml file(s).")

    # ---------------- ORGANISM-LEVEL ----------------
    print(f"\n{'#' * 60}\n# ORGANISM-LEVEL CONSERVATION ANALYSIS\n{'#' * 60}")
    org_summaries, org_occs = [], []
    for label in sorted(species_groups.keys()):
        summary_df, occ_df = run_conservation_analysis_for_group(
            label, species_groups[label], PLOT_DIR, tag_word="organism")
        if summary_df is not None:
            org_summaries.append(summary_df)
            org_occs.append(occ_df)

    if org_summaries:
        pd.concat(org_summaries, ignore_index=True).to_csv(OUT_SUMMARY_CSV, index=False)
        pd.concat(org_occs, ignore_index=True).to_csv(OUT_OCCURRENCE_CSV, index=False)
        print(f"\n[\u2713] ORGANISM-LEVEL: {len(org_summaries)} organism(s) analyzed -> "
              f"{OUT_SUMMARY_CSV}, {OUT_OCCURRENCE_CSV}")
        print(f"    Figures -> {PLOT_DIR}/<Organism>_motif_conservation_analysis.png")
    else:
        print("\n[!] ORGANISM-LEVEL: no organism produced usable motifs.")

    # ---------------- GENUS-LEVEL ----------------
    print(f"\n{'#' * 60}\n# GENUS-LEVEL CONSERVATION ANALYSIS\n{'#' * 60}")
    genus_summaries, genus_occs = [], []
    for label in sorted(genus_groups.keys()):
        summary_df, occ_df = run_conservation_analysis_for_group(
            label, genus_groups[label], PLOT_DIR_GENUS, tag_word="genus")
        if summary_df is not None:
            genus_summaries.append(summary_df)
            genus_occs.append(occ_df)

    if genus_summaries:
        pd.concat(genus_summaries, ignore_index=True).to_csv(OUT_SUMMARY_CSV_GENUS, index=False)
        pd.concat(genus_occs, ignore_index=True).to_csv(OUT_OCCURRENCE_CSV_GENUS, index=False)
        print(f"\n[\u2713] GENUS-LEVEL: {len(genus_summaries)} genus/genera analyzed -> "
              f"{OUT_SUMMARY_CSV_GENUS}, {OUT_OCCURRENCE_CSV_GENUS}")
        print(f"    Figures -> {PLOT_DIR_GENUS}/<Genus>_motif_conservation_analysis.png")
    else:
        print("\n[i] GENUS-LEVEL: no genus-level meme.xml files found (or no usable motifs) -- skipped.")


if __name__ == "__main__":
    main()
