#!/usr/bin/env python3
"""
STAGE 04d — Combined Protein Conservation Profile + Zoomed MSA Heatmap.

Generates TWO images in one run:

  1. protein_conservation_profile.png
     Wide, publication-style conservation curve across the full alignment
     (smoothed sliding-window %, matches the SARS-CoV-2-style reference
     figures). Fully vectorized -> fast and lightweight regardless of
     alignment size.

  2. msa_conservation_zoom.png
     A detailed, Clustal-colored per-residue heatmap, but only for a
     ZOOM window of the alignment (auto-selected as the most conserved
     region by default, or set manually via ZOOM_START / ZOOM_END below).
     This is what actually crashed before -- rendering a per-residue grid
     for a FULL alignment of hundreds of sequences x thousands of columns
     produces an unrenderable multi-hundred-megapixel image. Restricting
     it to a window (e.g. 150-300 columns) keeps it fast and readable.

Method notes:
  Per-column conservation = % of sequences matching the most common
  non-gap residue at that column (vectorized via NumPy boolean masks).
  Smoothing follows the standard sliding-window approach used by tools
  such as PlotSimilarity (GCG package) and AL2CO (Pei & Grishin, 2001,
  Bioinformatics 17(8):700-12).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO

ALIGNED_PATH = "outputs/family_aligned.fasta"
OUT_PROFILE_PATH = "outputs/protein_conservation_profile.png"
OUT_ZOOM_PATH = "outputs/msa_conservation_zoom.png"

GAP_CHAR = '-'
WINDOW_SIZE = 15            # smoothing window for the profile curve
CONSERVED_THRESHOLD = 80    # reference line, % conservation considered "highly conserved"

# Zoom panel settings. Leave both as None to auto-select the most conserved
# ZOOM_WIDTH-column window. Or set explicit 1-based alignment coordinates,
# e.g. ZOOM_START = 120, ZOOM_END = 320.
ZOOM_START = None
ZOOM_END = None
ZOOM_WIDTH = 200            # used only for auto-selection

# ClustalX Amino Acid Color Mapping (for the zoom heatmap panel)
CLUSTAL_COLORS = {
    'A': '#80a0f0', 'I': '#80a0f0', 'L': '#80a0f0', 'M': '#80a0f0', 'F': '#80a0f0', 'W': '#80a0f0', 'V': '#80a0f0',
    'R': '#f01505', 'K': '#f01505',
    'D': '#c048c0', 'E': '#c048c0',
    'N': '#15a015', 'Q': '#15a015', 'S': '#15a015', 'T': '#15a015',
    'C': '#f08080',
    'G': '#f08000',
    'P': '#c0c000',
    'H': '#15a0a0', 'Y': '#15a0a0',
    '-': '#ffffff'
}
DEFAULT_COLOR = '#e0e0e0'


def _hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))


def alignment_to_byte_array(alignment):
    """Converts alignment to a (num_seqs x align_len) uint8 array of ASCII codes."""
    seqs = [str(rec.seq).upper() for rec in alignment]
    align_len = len(seqs[0])
    for i, s in enumerate(seqs):
        if len(s) != align_len:
            sys.exit(
                f"[!] Error: sequence '{alignment[i].id}' has length {len(s)}, "
                f"expected {align_len}. Alignment columns must be equal length."
            )
    seq_bytes = np.array(
        [np.frombuffer(s.encode('ascii'), dtype=np.uint8) for s in seqs],
        dtype=np.uint8
    )
    return seq_bytes


def compute_raw_conservation(seq_bytes, num_seqs):
    """Vectorized per-column conservation score (%)."""
    gap_byte = ord(GAP_CHAR)
    unique_chars = np.unique(seq_bytes)
    align_len = seq_bytes.shape[1]
    max_freq = np.zeros(align_len, dtype=np.int64)

    for uc in unique_chars:
        if uc == gap_byte:
            continue
        counts = (seq_bytes == uc).sum(axis=0)
        max_freq = np.maximum(max_freq, counts)

    return (max_freq / num_seqs) * 100


def smooth_conservation(scores, window_size):
    """Sliding-window running average (PlotSimilarity/AL2CO convention)."""
    if window_size <= 1:
        return scores
    kernel = np.ones(window_size) / window_size
    return np.convolve(scores, kernel, mode='same')


def compute_gap_frequency(seq_bytes, num_seqs):
    gap_byte = ord(GAP_CHAR)
    gap_counts = (seq_bytes == gap_byte).sum(axis=0)
    return (gap_counts / num_seqs) * 100


def auto_select_zoom_window(smoothed_scores, width):
    """Finds the WIDTH-column window with the highest average smoothed conservation."""
    align_len = len(smoothed_scores)
    width = min(width, align_len)
    cumsum = np.cumsum(np.insert(smoothed_scores, 0, 0))
    window_sums = cumsum[width:] - cumsum[:-width]
    best_start = int(np.argmax(window_sums))
    return best_start, best_start + width  # 0-based [start, end)


def plot_conservation_profile(x_positions, smoothed_scores, gap_freq, num_seqs, align_len,
                               zoom_range, out_path):
    fig, ax = plt.subplots(figsize=(14, 5), dpi=200)

    ax.plot(x_positions, smoothed_scores, color='#1f4fd8', linewidth=1.3,
            label=f'Conservation (smoothed, window={WINDOW_SIZE})')
    ax.fill_between(x_positions, smoothed_scores, color='#1f4fd8', alpha=0.25)
    ax.axhline(CONSERVED_THRESHOLD, color='firebrick', linestyle='--',
               linewidth=1, alpha=0.7, label=f'{CONSERVED_THRESHOLD}% conservation threshold')

    zoom_start_1based = zoom_range[0] + 1
    zoom_end_1based = zoom_range[1]
    ax.axvspan(zoom_start_1based, zoom_end_1based, color='orange', alpha=0.15,
               label=f'Zoom region ({zoom_start_1based}-{zoom_end_1based})')

    ax.set_xlim(1, align_len)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Alignment Position (aa)", fontsize=11)
    ax.set_ylabel("Conservation (%)", fontsize=11)
    ax.set_title(
        f"Protein Conservation Profile ({num_seqs} sequences, {align_len} aligned positions)",
        fontsize=13, pad=12
    )
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
    ax.grid(True, linestyle='--', alpha=0.4)

    ax2 = ax.twinx()
    ax2.plot(x_positions, gap_freq, color='gray', linewidth=0.6, alpha=0.35)
    ax2.set_ylabel("Gap frequency (%)", fontsize=9, color='gray')
    ax2.set_ylim(0, 105)
    ax2.tick_params(axis='y', labelsize=7, colors='gray')

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"[✓] Conservation profile plot saved -> {out_path}")


def plot_zoom_heatmap(alignment, seq_bytes, zoom_range, out_path):
    """Detailed Clustal-colored per-residue heatmap, restricted to zoom_range columns."""
    start, end = zoom_range
    sub_bytes = seq_bytes[:, start:end]
    num_seqs, width = sub_bytes.shape

    palette = np.full((256, 3), _hex_to_rgb(DEFAULT_COLOR), dtype=np.uint8)
    for aa, hex_color in CLUSTAL_COLORS.items():
        palette[ord(aa)] = _hex_to_rgb(hex_color)
    img = palette[sub_bytes]

    fig_height = max(6, num_seqs * 0.16 + 2)
    fig, ax = plt.subplots(figsize=(max(10, width * 0.12), fig_height), dpi=200)

    ax.imshow(
        img,
        extent=(start + 0.5, end + 0.5, 0.5, num_seqs + 0.5),
        origin='upper',
        aspect='auto',
        interpolation='none'
    )

    if width <= 250:
        for row_idx, record in enumerate(alignment):
            seq = str(record.seq).upper()[start:end]
            y = num_seqs - row_idx
            for col_idx, aa in enumerate(seq):
                ax.text(start + col_idx + 1, y, aa, ha='center', va='center',
                        fontsize=6, fontweight='bold')

    ax.set_xlim(start + 0.5, end + 0.5)
    ax.set_ylim(0.5, num_seqs + 0.5)
    ax.set_yticks(range(1, num_seqs + 1))
    ax.set_yticklabels([rec.id for rec in reversed(alignment)], fontsize=7)
    ax.set_xlabel("Alignment Position (aa)", fontsize=10)
    ax.set_title(f"Zoomed MSA View (positions {start + 1}-{end})", fontsize=12, pad=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"[✓] Zoomed MSA heatmap saved -> {out_path}")


def main():
    if not (os.path.exists(ALIGNED_PATH) and os.path.getsize(ALIGNED_PATH) > 0):
        sys.exit(f"[!] Error: Alignment file '{ALIGNED_PATH}' not found or empty.")

    alignment = AlignIO.read(ALIGNED_PATH, "fasta")
    num_seqs = len(alignment)
    seq_bytes = alignment_to_byte_array(alignment)
    align_len = seq_bytes.shape[1]
    print(f"[*] Alignment loaded: {num_seqs} sequences x {align_len} columns")

    print("[*] Computing per-column conservation (vectorized)...")
    raw_scores = compute_raw_conservation(seq_bytes, num_seqs)

    print(f"[*] Smoothing with sliding window (size={WINDOW_SIZE})...")
    smoothed_scores = smooth_conservation(raw_scores, WINDOW_SIZE)

    print("[*] Computing gap frequency track...")
    gap_freq = compute_gap_frequency(seq_bytes, num_seqs)

    x_positions = np.arange(1, align_len + 1)

    if ZOOM_START is not None and ZOOM_END is not None:
        zoom_range = (ZOOM_START - 1, ZOOM_END)  # convert to 0-based [start, end)
        print(f"[*] Using manual zoom window: {ZOOM_START}-{ZOOM_END}")
    else:
        zoom_range = auto_select_zoom_window(smoothed_scores, ZOOM_WIDTH)
        print(f"[*] Auto-selected most conserved zoom window: "
              f"{zoom_range[0] + 1}-{zoom_range[1]}")

    plot_conservation_profile(x_positions, smoothed_scores, gap_freq, num_seqs, align_len,
                               zoom_range, OUT_PROFILE_PATH)

    plot_zoom_heatmap(alignment, seq_bytes, zoom_range, OUT_ZOOM_PATH)

    print("\n[✓] Done. Two images generated:")
    print(f"    1) {OUT_PROFILE_PATH}")
    print(f"    2) {OUT_ZOOM_PATH}")


if __name__ == "__main__":
    main()