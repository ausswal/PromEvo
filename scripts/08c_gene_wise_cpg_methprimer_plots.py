#!/usr/bin/env python3
"""
STAGE 08c — Gene-wise CpG Island Visualization (MethPrimer-style).

Direct Python translation of CpG_MethPrimer_COUG_code.R, wired to this
pipeline's own data instead of a hardcoded local Windows path. Produces ONE
three-panel PNG per individual promoter sequence (Obs/Exp CpG ratio, GC%
with shaded putative islands, and a schematic islands/CpG-dinucleotide
track) -- exactly the classic MethPrimer figure, not the adaptive-threshold
summary analysis in 08_CpG_Island_adjusted_Finder_for_all.py /
08_new_organism_wise_cpg_island_finder.py.

THIS IS A DIFFERENT ALGORITHM FROM STAGE 08, ON PURPOSE
---------------------------------------------------------------------------
Stage 08's adaptive threshold picks a GC%/Obs-Exp cutoff FROM the data
(per-organism or per-genus global GC%, Gardiner-Garden vs AT-rich
fallback). This script instead reproduces the classic FIXED Gardiner &
Garden criteria the R script used (GC% >= 50, Obs/Exp CpG >= 0.60,
sliding 100 bp window), matching what MethPrimer itself reports -- these
two scripts can and will disagree on which regions count as "islands" for
the same sequence, particularly for very GC-rich or GC-poor organisms.
That's expected: they answer different questions (data-adaptive regional
comparison vs. the canonical fixed-threshold definition), not a bug in
either one. Do not treat this script's output as "the corrected version of"
Stage 08's, or vice versa.

METHOD LABEL FOR YOUR METHODS SECTION: this is the "CLASSIC" CpG-island method
(fixed GC >= 50%, Obs/Exp >= 0.6, 100-bp window). Stage 08 is the "ADAPTIVE"
method (200-bp window, GC cutoff derived from the group's global GC%).

COORDINATES / FULL-LENGTH PROFILE (changed after review)
---------------------------------------------------------------------------
* The x-axis is TSS-RELATIVE like the rest of the pipeline: each promoter's
  3' end is the TSS, so the last base is -1 and a 1000-bp promoter spans
  -1000 ... -1.
* Each window's value is plotted at the window CENTRE, and the curves are
  extended flat to both ends (same "fill=extend" idea as rollmean_extend), so
  the whole promoter is drawn -- previously only window starts 1..901 were
  shown, dropping the 99 bp closest to the TSS and clipping CpG ticks there.
* An island is shaded over every base covered by its qualifying windows.

WHAT CHANGED FROM THE R VERSION (beyond the language)
---------------------------------------------------------------------------
1. Input/output paths: the R script had one hardcoded fasta_file/outdir for
   a single local promoter set. This reads outputs/promoter_sequences.fasta
   (the same master promoter FASTA every other stage in this pipeline
   uses) and organizes output PNGs by organism, using outputs/
   family_labels.csv's "organism" column and the same sanitize_organism_
   name() convention as 06b/06c/08-organism/23_comparative_matrix.py, so a
   given organism's folder name matches everywhere in this pipeline.
2. Panel 3's x-axis ticks: the R script hardcoded `bp <- seq(0,900,200)`
   regardless of the actual sequence length -- fine for a dataset of
   uniform ~1000 bp promoters, but wrong for anything shorter/longer (and
   this pipeline's PROMOTER_UPSTREAM_BP is user-configurable at fetch
   time). This version computes ticks from the actual sequence length,
   same as Panels 1 and 2 already did in the R original.
3. The nested-loop GC%/Obs-Exp computation (one Python-level loop per
   window position, itself doing 3 substring scans) is replaced with a
   vectorized sliding-window sum via cumulative sums -- verified to
   produce bit-identical results to the naive translation on test data,
   just fast enough to run across an entire promoter set instead of one
   file at a time.
4. zoo::rollmean(..., fill="extend") (centered window=5 moving average,
   edge positions filled by repeating the nearest full-window value) is
   reproduced exactly via a centered convolution + edge-replication -- see
   rollmean_extend() below.
5. R's `str_count(s, "CG")` / `gregexpr("CG", seq)` and Python's
   str.count("CG") / re.finditer("CG", seq) are equivalent here: "CG" as a
   2-character pattern cannot self-overlap (unlike e.g. "AA"), so
   overlapping vs. non-overlapping counting gives identical results in
   both languages -- confirmed, not assumed.

Requirements:
    pip install biopython numpy matplotlib pandas
"""

import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless-safe; this script only ever saves PNGs
import matplotlib.pyplot as plt
from Bio import SeqIO

from pipeline_config import check_promoter_lengths

# ============================== CONFIGURATION ==============================

PROMOTER_FASTA_PATH = "outputs/promoter_sequences.fasta"   # same source every other stage uses
LABELS_PATH = "outputs/family_labels.csv"                   # confirmed "organism" column (Stage 01)

OUT_ROOT = "outputs/CpG_Island_GeneWise"

WINDOW = 100          # sliding window size, bp (Gardiner & Garden)
ROLL_K = 5             # centered smoothing window for GC%/Obs-Exp curves
GC_CUTOFF = 50.0       # % GC threshold
OE_CUTOFF = 0.60       # Observed/Expected CpG threshold

# Matches the R script's png(width=2800, height=3200, res=600) exactly.
FIG_DPI = 600
FIG_WIDTH_IN = 2800 / FIG_DPI
FIG_HEIGHT_IN = 3200 / FIG_DPI

GRIDLINE_COLOR = "#D8FFFF"
ISLAND_FILL = "#A8DFF0"
ISLAND_BORDER = "#3F88C5"
CPG_TICK_COLOR = "#FF8C00"

# ===========================================================================


def sanitize_organism_name(name):
    """Same alnum-only sanitizer used throughout the pipeline (01, 06b/06c,
    08-organism, 23_comparative_matrix.py) so folder names line up
    everywhere."""
    if not name or str(name).strip() == "" or str(name).lower() == "nan":
        return "Unknown_organism"
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_")


def load_organism_lookup(labels_path):
    """Returns {seq_id: organism}. Missing file / missing column is a soft
    failure -- every sequence just falls into 'Unknown_organism' rather
    than crashing, since gene-wise plots are still meaningful without
    organism grouping."""
    if not os.path.exists(labels_path):
        print(f"[!] {labels_path} not found -- all gene plots will be grouped "
              f"under 'Unknown_organism'.")
        return {}
    labels_df = pd.read_csv(labels_path)
    if "organism" not in labels_df.columns or "seq_id" not in labels_df.columns:
        print(f"[!] {labels_path} is missing 'seq_id'/'organism' columns -- "
              f"all gene plots will be grouped under 'Unknown_organism'.")
        return {}
    return dict(zip(labels_df["seq_id"], labels_df["organism"]))


def sliding_window_sum(bool_arr, window):
    """Sum of a boolean/0-1 array over every window-length span, via
    cumulative sums. Returns an array of length len(bool_arr) - window + 1."""
    csum = np.cumsum(np.insert(bool_arr.astype(np.int64), 0, 0))
    return csum[window:] - csum[:-window]


def compute_gc_oe_profile(seq, window=WINDOW):
    """
    Vectorized equivalent of the R script's per-window loop:
        C, G, CG counted in each window-length substring
        GC[i]  = (C+G)/window*100
        OE[i]  = 0 if C==0 or G==0 else (CG*window)/(C*G)
    Returns (GC, OE), each length len(seq) - window + 1.
    """
    n = len(seq)
    arr = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)

    C_win = sliding_window_sum(arr == ord("C"), window)
    G_win = sliding_window_sum(arr == ord("G"), window)

    # CG dinucleotide flags, one per possible start position (length n-1).
    # A window starting at i covers CG-start-positions i .. i+window-2, i.e.
    # a sliding sum of cg_flags with window length (window-1).
    cg_flags = (arr[:-1] == ord("C")) & (arr[1:] == ord("G"))
    cg_win = sliding_window_sum(cg_flags, window - 1)

    GC = (C_win + G_win) / window * 100.0
    with np.errstate(divide="ignore", invalid="ignore"):
        OE = np.where((C_win == 0) | (G_win == 0), 0.0, (cg_win * window) / (C_win * G_win))

    return GC, OE


def rollmean_extend(x, k=ROLL_K):
    """
    Reproduces zoo::rollmean(x, k, fill="extend"): a centered moving
    average that, for k=5, loses the outer 2 points on each side to
    edge effects, then fills those edge positions by repeating the
    nearest fully-computed value (rather than leaving them NA).
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n == 0:
        return x
    if n < k:
        return np.full(n, x.mean())

    kernel = np.ones(k) / k
    valid = np.convolve(x, kernel, mode="valid")   # length n - k + 1
    half = (k - 1) // 2

    out = np.empty(n)
    out[half:n - half] = valid
    out[:half] = valid[0]
    out[n - half:] = valid[-1]
    return out


def find_runs(bool_arr):
    """Run-length encoding of a boolean array. Returns a list of
    (start_idx, end_idx, value) with 0-based, inclusive indices --
    equivalent to R's rle() + cumsum-derived starts/ends."""
    if len(bool_arr) == 0:
        return []
    runs = []
    start = 0
    current = bool_arr[0]
    for i in range(1, len(bool_arr)):
        if bool_arr[i] != current:
            runs.append((start, i - 1, current))
            start = i
            current = bool_arr[i]
    runs.append((start, len(bool_arr) - 1, current))
    return runs


def plot_gene_cpg(seq_id, seq, out_path, strain_title):
    """One three-panel MethPrimer-style figure for a single sequence, drawn
    over the FULL promoter on a TSS-relative axis (last base = -1)."""

    n = len(seq)
    if n < WINDOW:
        print(f"    [!] {seq_id}: length {n} bp < WINDOW={WINDOW} bp -- skipped.")
        return False

    GC_raw, OE_raw = compute_gc_oe_profile(seq, WINDOW)
    GC = rollmean_extend(GC_raw, ROLL_K)
    OE = rollmean_extend(OE_raw, ROLL_K)
    n_win = len(GC)                                  # == n - WINDOW + 1

    island_mask = (GC >= GC_CUTOFF) & (OE >= OE_CUTOFF)
    island_runs = [(st, en) for st, en, val in find_runs(island_mask) if val]

    # --- coordinates -------------------------------------------------------
    # Base j (0-based) sits at TSS-relative x = j - n  (last base -> -1).
    # Window starting at s covers bases s .. s+WINDOW-1; its centre is
    # s + (WINDOW-1)/2, so it is drawn at that centre.
    centres = np.arange(n_win) + (WINDOW - 1) / 2.0 - n
    x_left, x_right = -float(n), -1.0
    base_x = np.arange(n) - n                        # one x per base

    def extend(y):
        """Window-centre curve -> curve spanning the whole promoter (flat
        extension past the first/last window centre)."""
        return (np.concatenate(([x_left], centres, [x_right])),
                np.concatenate(([y[0]], y, [y[-1]])))

    xe, GCe = extend(GC)
    _, OEe = extend(OE)

    # Bases covered by at least one qualifying window (for shading).
    base_island = np.zeros(n, dtype=bool)
    for st, en in island_runs:
        base_island[st:en + WINDOW] = True           # windows st..en -> bases st..en+WINDOW-1

    # CpG dinucleotide ticks: centre of each "CG" (bases j, j+1), whole sequence.
    cpg_x = np.array([m.start() + 0.5 - n for m in re.finditer("CG", seq)])

    fig, axes = plt.subplots(
        3, 1, figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        gridspec_kw={"hspace": 1.15},
    )
    fig.subplots_adjust(top=0.85, bottom=0.09, left=0.16, right=0.96)
    fig.suptitle(strain_title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(0.5, 0.925,
             f"Classic CpG island criteria: {WINDOW}-bp window, GC \u2265 {GC_CUTOFF:g}%, "
             f"Obs/Exp \u2265 {OE_CUTOFF:g}",
             ha="center", va="center", fontsize=8.5, color="#444444")

    # Round TSS-relative ticks (-1000, -800, ... -200) plus -1 at the TSS end;
    # the left edge (-n) is the axis limit itself, so clipped promoters stay uncluttered.
    tick_positions = np.append(np.arange(-(n // 200) * 200, 0, 200), -1)
    tick_labels = [str(int(t)) for t in tick_positions]
    grid_positions = np.arange(-n, 0, 10)
    xlabel = "Position relative to TSS (bp)"

    def style_axis(ax):
        for gx in grid_positions:
            ax.axvline(gx, color=GRIDLINE_COLOR, linewidth=0.5, zorder=0)
        ax.set_xlim(x_left, x_right)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=9)
        ax.set_xlabel(xlabel, fontsize=11, fontweight="bold")

    # ---------------- PANEL 1: Observed vs Expected ----------------
    ax1 = axes[0]
    style_axis(ax1)
    ax1.axhline(OE_CUTOFF, color="darkorange", linestyle="--", linewidth=2, zorder=2)
    ax1.plot(xe, OEe, color="black", linewidth=2.5, zorder=3)
    ax1.set_ylim(0, max(1.5, float(np.max(OE))))
    ax1.set_ylabel("Obs/Exp", fontsize=11, fontweight="bold")
    ax1.set_title("Observed vs Expected", fontsize=11, fontweight="bold", pad=8)

    # ---------------- PANEL 2: GC Percentage ----------------
    ax2 = axes[1]
    style_axis(ax2)
    ax2.axhline(GC_CUTOFF, color="red", linewidth=2, linestyle="--", zorder=2)
    gc_base = np.interp(base_x, xe, GCe)
    ax2.fill_between(base_x, gc_base, 0, where=base_island,
                     color=ISLAND_FILL, linewidth=0, zorder=1, interpolate=True)
    ax2.plot(xe, GCe, color="black", linewidth=2.5, zorder=3)
    ax2.set_ylim(0, 100)
    ax2.set_ylabel("GC Percentage", fontsize=11, fontweight="bold")
    ax2.set_title("GC Percentage", fontsize=12, fontweight="bold", pad=6)

    # ---------------- PANEL 3: Putative Islands ----------------
    ax3 = axes[2]
    style_axis(ax3)
    ax3.plot([x_left, x_right], [0, 0], color="red", linewidth=5, zorder=2)
    for st, en in island_runs:
        x0 = st - n - 0.5                            # first covered base (left edge)
        x1 = (en + WINDOW - 1) - n + 0.5             # last covered base (right edge)
        ax3.add_patch(plt.Rectangle(
            (x0, 0), x1 - x0, 1,
            facecolor=ISLAND_FILL, edgecolor=ISLAND_BORDER, linewidth=1, zorder=1,
        ))
    if len(cpg_x):
        ax3.vlines(cpg_x, -0.08, 0.08, color=CPG_TICK_COLOR, linewidth=1, zorder=3)
    ax3.set_ylim(0, 1)
    ax3.set_ylabel("Threshold", fontsize=11, fontweight="bold")
    ax3.set_title("Putative Islands", fontsize=12, fontweight="bold", pad=6)

    fig.savefig(out_path, dpi=FIG_DPI)
    plt.close(fig)
    return True


def main():
    if not os.path.exists(PROMOTER_FASTA_PATH):
        sys.exit(f"[!] {PROMOTER_FASTA_PATH} not found -- run the promoter-fetch stage first.")

    os.makedirs(OUT_ROOT, exist_ok=True)

    organism_lookup = load_organism_lookup(LABELS_PATH)

    records = list(SeqIO.parse(PROMOTER_FASTA_PATH, "fasta"))
    if not records:
        sys.exit(f"[!] No sequences parsed from {PROMOTER_FASTA_PATH}.")

    print("=" * 65)
    print("  GENE-WISE CpG ISLAND VISUALIZATION (MethPrimer-style)  ")
    print("=" * 65)
    print(f"[*] {len(records)} promoter sequence(s) to plot")
    check_promoter_lengths([len(r.seq) for r in records], stage="08c")
    print(f"[*] WINDOW={WINDOW} bp, GC_CUTOFF={GC_CUTOFF}%, OE_CUTOFF={OE_CUTOFF}")
    print()

    n_plotted, n_skipped = 0, 0

    for rec in records:
        seq_id = rec.id
        seq = str(rec.seq).upper()

        organism = organism_lookup.get(seq_id)
        org_label = sanitize_organism_name(organism)

        gene = re.sub(r"[^A-Za-z0-9_]", "_", seq_id)
        parts = gene.split("_")
        strain_title = f"{parts[0]}_{parts[1]}" if len(parts) >= 2 else gene

        org_dir = os.path.join(OUT_ROOT, org_label)
        os.makedirs(org_dir, exist_ok=True)
        out_path = os.path.join(org_dir, f"{gene}.png")

        ok = plot_gene_cpg(seq_id, seq, out_path, strain_title)
        n_plotted += int(ok)
        n_skipped += int(not ok)

    print()
    print(f"[\u2713] Gene-wise CpG island plots complete: {n_plotted} plotted, "
          f"{n_skipped} skipped (sequence shorter than WINDOW={WINDOW} bp).")
    print(f"    Outputs under: {OUT_ROOT}/<Organism>/<gene>.png")


if __name__ == "__main__":
    main()
