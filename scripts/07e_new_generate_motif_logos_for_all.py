#!/usr/bin/env python3
"""
STAGE 16 / 7e — Enriched Motif Logo Generator (Single Frame Visual).

Features:
1. Parses MEME XML output (meme.xml) to extract Position Probability Matrices (PPMs) 
   and calculates Information Content (Bits per position).
2. Uses Matplotlib PathPatch vector affine transformations to stretch nucleotide glyphs 
   proportional to bit height (0.0 to 2.0 bits).
3. Combines all enriched motif logos into a clean multi-panel stacked frame.
4. Exports a high-resolution, publication-ready figure.
"""

import os
import re
import sys
import glob
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import matplotlib.transforms as transforms

# Try importing logomaker if installed
try:
    import logomaker
    HAS_LOGOMAKER = True
except ImportError:
    HAS_LOGOMAKER = False

# ============================== CONFIGURATION ==============================
MEME_OUT_DIR = "outputs/meme_output"                # species-level (07c)
MEME_OUT_DIR_GENUS = "outputs/meme_output_genus"     # genus-level (07c)
LEGACY_MEME_XML = "outputs/meme_output/meme.xml"     # pre-split fallback

OUT_PLOT_DIR_ORG = "outputs/motif_logos_by_organism"
OUT_PLOT_DIR_GENUS = "outputs/motif_logos_by_genus"
COLOR_SCHEME = {'A': '#2CA02C', 'C': '#1F77B4', 'G': '#FF7F0E', 'T': '#D62728'}
# ===========================================================================


def sanitize_filename(name):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip())
    return safe.strip("_") or "Unknown"


def discover_meme_xml_groups(meme_out_dir):
    """Finds every <meme_out_dir>/<Label>/meme.xml and returns {label: path}."""
    groups = {}
    for xml_path in sorted(glob.glob(os.path.join(meme_out_dir, "*", "meme.xml"))):
        label = os.path.basename(os.path.dirname(xml_path))
        groups[label] = xml_path
    return groups


def parse_meme_matrices(xml_path):
    """Parses motif PPMs (Position Probability Matrices) from meme.xml."""
    if not os.path.exists(xml_path):
        print(f"[!] Error: Could not find MEME XML file at: {xml_path}")
        return None

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        motifs_dict = {}
        for motif in root.findall(".//motifs/motif"):
            m_id = motif.get("id")
            m_name = motif.get("name")
            e_value = float(motif.get("e_value"))
            alt = motif.get("alt") or m_name

            # Extract probability matrix
            matrix_data = []
            for pos in motif.findall(".//probabilities/alphabet_matrix/alphabet_array"):
                probs = [float(child.text) for child in pos.findall("value")]
                matrix_data.append(probs)

            # Standard ACGT order in MEME XML
            df_ppm = pd.DataFrame(matrix_data, columns=["A", "C", "G", "T"])
            
            # Convert PPM to Information Content (PWM in bits)
            # Bits = 2.0 + sum(p * log2(p))
            entropy = -np.sum(np.nan_to_num(df_ppm * np.log2(df_ppm + 1e-9)), axis=1)
            info_content = 2.0 - entropy
            df_pwm = df_ppm.multiply(info_content, axis=0)

            motifs_dict[m_id] = {
                "name": m_name,
                "alt": alt,
                "e_value": e_value,
                "ppm": df_ppm,
                "pwm": df_pwm
            }

        return motifs_dict
    except Exception as e:
        print(f"Warning: Failed to parse {xml_path}: {e}")
        return None


def draw_vector_glyph(ax, letter, x, y, width, height, color):
    """Draws a letter vector glyph transformed precisely into [x, x+width] x [y, y+height]."""
    if height <= 0.005:
        return
        
    prop = FontProperties(family='sans-serif', weight='bold')
    tp = TextPath((0, 0), letter, size=1, prop=prop)
    bbox = tp.get_extents()
    
    scale_x = width / (bbox.width if bbox.width > 0 else 1)
    scale_y = height / (bbox.height if bbox.height > 0 else 1)
    
    trans = transforms.Affine2D() \
        .translate(-bbox.x0, -bbox.y0) \
        .scale(scale_x, scale_y) \
        .translate(x, y) + ax.transData
        
    patch = PathPatch(tp, facecolor=color, edgecolor='black', linewidth=0.3, transform=trans)
    ax.add_patch(patch)


def plot_custom_vector_logo(ax, df_pwm):
    """Custom vector motif rendering engine using affine-transformed glyphs."""
    ax.set_ylim(0, 2.0)
    ax.set_xlim(0.5, len(df_pwm) + 0.5)
    
    for pos_idx, row in df_pwm.iterrows():
        x_pos = pos_idx + 1
        sorted_bases = row.sort_values()
        y_bottom = 0.0
        
        for base, height in sorted_bases.items():
            if height > 0.005:
                draw_vector_glyph(
                    ax=ax, 
                    letter=base, 
                    x=x_pos - 0.35, 
                    y=y_bottom, 
                    width=0.70, 
                    height=height, 
                    color=COLOR_SCHEME[base]
                )
                y_bottom += height

    ax.set_xticks(range(1, len(df_pwm) + 1))
    ax.set_xticklabels(range(1, len(df_pwm) + 1), fontsize=8.5, fontweight='bold')


def generate_logo_frame_for_group(label, xml_path, out_dir, tag_word):
    """Builds the multi-panel stacked logo figure for ONE group's (one
    organism's or one genus's) meme.xml. Returns True on success."""
    print(f"\n{'-' * 60}\n{tag_word.capitalize()}: {label}\n{'-' * 60}")

    motifs_data = parse_meme_matrices(xml_path)
    if not motifs_data:
        print(f"    [!] No motifs parsed from {xml_path} — skipping.")
        return False

    num_motifs = len(motifs_data)
    fig, axes = plt.subplots(
        nrows=num_motifs, ncols=1, figsize=(11, 2.2 * num_motifs), dpi=200, sharex=False
    )
    if num_motifs == 1:
        axes = [axes]

    for idx, (m_id, m_info) in enumerate(motifs_data.items()):
        ax = axes[idx]
        df_pwm = m_info["pwm"]

        if HAS_LOGOMAKER:
            logo = logomaker.Logo(
                df_pwm, ax=ax, color_scheme=COLOR_SCHEME, baseline_width=1, show_spines=True
            )
            logo.style_spines(visible=True, linewidth=1, color='black')
            logo.style_xticks(fmt='%d', size=9, fontweight='bold')
            ax.set_ylabel("Bits", fontsize=9.5, fontweight='bold')
            ax.set_ylim(0, 2.0)
        else:
            plot_custom_vector_logo(ax, df_pwm)
            ax.set_ylabel("Bits", fontsize=9.5, fontweight='bold')

        title_text = (f"Motif {idx + 1} ({m_info['alt']})  |  E-value: {m_info['e_value']:.2e}  |  "
                      f"Length: {len(df_pwm)} bp")
        ax.set_title(title_text, fontsize=10.5, fontweight='bold', loc='left', pad=4, color='#1A252C')
        ax.set_xlabel("Position (bp)", fontsize=9, fontweight='bold')

    fig_height = 2.2 * num_motifs
    plt.suptitle(f"Enriched Motif Sequence Logos — {tag_word.capitalize()}: {label}",
                 fontsize=12.5, fontweight='bold', y=0.99)
    plt.tight_layout()
    # Reserve a fixed ~0.55in top margin for the suptitle regardless of how
    # many motifs (and therefore how tall the figure) there are -- a fixed
    # top FRACTION (e.g. 0.93) leaves plenty of room for many-motif figures
    # but collides with the top panel's own title when there's only 1-2
    # motifs, since the absolute margin shrinks with the fraction.
    top_margin_in = 0.55
    plt.subplots_adjust(top=max(0.75, 1 - top_margin_in / fig_height))

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sanitize_filename(label)}_logos.png")
    plt.savefig(out_path, dpi=200, bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)

    print(f"    [✓] {num_motifs} motif(s) -> {out_path}")
    return True


def main():
    os.makedirs("outputs", exist_ok=True)

    species_groups = discover_meme_xml_groups(MEME_OUT_DIR)
    genus_groups = discover_meme_xml_groups(MEME_OUT_DIR_GENUS)

    if not species_groups and not genus_groups:
        if os.path.exists(LEGACY_MEME_XML):
            print(f"[!] No per-organism/per-genus meme.xml files found under {MEME_OUT_DIR}*/ — "
                  f"falling back to the single legacy file at {LEGACY_MEME_XML}.")
            species_groups = {"All_organisms": LEGACY_MEME_XML}
        else:
            sys.exit(f"[!] No meme.xml files found under '{MEME_OUT_DIR}/<Organism>/', "
                      f"'{MEME_OUT_DIR_GENUS}/<Genus>/', or '{LEGACY_MEME_XML}'. "
                      f"Run 07c_Motif_discovery_for_all.py first.")

    print(f"[+] Found {len(species_groups)} organism-level and {len(genus_groups)} "
          f"genus-level meme.xml file(s).")

    print(f"\n{'#' * 60}\n# SPECIES-LEVEL LOGOS\n{'#' * 60}")
    n_species_ok = 0
    for label in sorted(species_groups.keys()):
        if generate_logo_frame_for_group(label, species_groups[label], OUT_PLOT_DIR_ORG, "organism"):
            n_species_ok += 1

    print(f"\n{'#' * 60}\n# GENUS-LEVEL LOGOS\n{'#' * 60}")
    n_genus_ok = 0
    for label in sorted(genus_groups.keys()):
        if generate_logo_frame_for_group(label, genus_groups[label], OUT_PLOT_DIR_GENUS, "genus"):
            n_genus_ok += 1

    print(f"\n[✓] Generated {n_species_ok}/{len(species_groups)} organism-level and "
          f"{n_genus_ok}/{len(genus_groups)} genus-level logo figure(s).")
    print(f"    Species figures -> {OUT_PLOT_DIR_ORG}/<Organism>_logos.png")
    print(f"    Genus figures   -> {OUT_PLOT_DIR_GENUS}/<Genus>_logos.png")


if __name__ == "__main__":
    main()