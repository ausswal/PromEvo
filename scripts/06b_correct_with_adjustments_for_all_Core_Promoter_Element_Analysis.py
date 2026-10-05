#!/usr/bin/env python3
"""
STAGE 06b — Core Promoter Element Analysis for UPSTREAM-ONLY promoters

Reads promoter sequences from outputs/promoter_sequences.fasta (the 1000 bp
immediately upstream of the annotated TSS / gene start; last base = -1; no
downstream flank is analysed anywhere in this pipeline):
1. Scans promoters for the core motifs (TATA, Inr, DPE, BREu, MTE, CAAT) as
   IUPAC sequence patterns on BOTH forward (+) and reverse (-) strands.
2. Maps each hit's START to a TSS-relative position (rel_pos = start - seq_len,
   so the last base is -1; end-anchored, so promoters clipped shorter than
   1000 bp at a contig edge are still aligned on the TSS).
3. Bins hits into UPSTREAM-ONLY positional windows:
       Distal    (-1000 to -151)
       Proximal  (-150  to -31)
       Core      (-30   to -1)
   (v2 FIX: the previous "-30 to +20" and "+21 to +200" windows can never be
   populated by an upstream-only promoter -- the largest possible rel_pos is
   -1 -- so the downstream bin was always 0 and the core bin was really only
   -30..-1. The windows now say exactly what they measure.)
4. Exports a summary CSV and a granular individual-hits CSV.
5. Evaluates signal-to-noise (SNR) against composition-shuffled controls.
6. Draws the 4-panel figure (plot_core_promoter_figure -- shared with
   06b-organism so the two cannot drift apart).

NOTE ON DPE / MTE / Inr: canonically these elements sit at or just
downstream of the TSS (Inr ~ -2..+4, MTE ~ +18..+29, DPE ~ +28..+33). In an
upstream-only region, hits are therefore sequence-pattern matches in the
upstream region, interpreted via SNR against the shuffled background -- not
positioned core elements. The figure carries a one-line footnote saying so.
"""

import os
import csv
import random
import re
import numpy as np
import matplotlib.pyplot as plt

INPUT_FASTA = "outputs/promoter_sequences.fasta"
OUT_PLOT = "outputs/plasmodium_promoter_conservation_scan.png"
OUT_CSV_SUMMARY = "outputs/core_promoter_summary_metrics.csv"
OUT_CSV_HITS = "outputs/core_promoter_individual_hits.csv"

# Stage 01 coordinates: promoters are the 1000 bp upstream of the TSS, no downstream flank
UPSTREAM_LEN = 1000
DOWNSTREAM_LEN = 0
FULL_LEN = UPSTREAM_LEN + DOWNSTREAM_LEN  # 1000 bp total

random.seed(42)
np.random.seed(42)

# Upstream-only positional windows (TSS-relative; last base = -1). Membership is by the
# motif's START position. (key, low, high, display label)
WINDOWS = [
    ("distal",   -10**9,       -151, f"Distal (-{UPSTREAM_LEN} to -151)"),
    ("proximal", -150,         -31,  "Proximal (-150 to -31)"),
    ("core",     -30,          -1,   "Core (-30 to -1)"),
]
UPSTREAM_ONLY_NOTE = ("Upstream-only promoters (-%d to -1): DPE/MTE/Inr canonically lie at/after the TSS, "
                      "so their hits here are sequence-pattern matches judged against the shuffled "
                      "background (SNR), not positioned core elements." % UPSTREAM_LEN)


def window_of(rel_pos):
    for key, lo, hi, _ in WINDOWS:
        if lo <= rel_pos <= hi:
            return key
    return None

# IUPAC ambiguity helper to compile regex patterns
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': '[AG]', 'Y': '[CT]', 
    'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 
    'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

def reverse_complement(seq):
    """Generates reverse complement of a nucleotide sequence."""
    return seq.translate(COMPLEMENT)[::-1]

def parse_iupac(pat):
    return re.compile("".join(IUPAC.get(b, b) for b in pat.upper()))

# Core Eukaryotic & Plasmodium Promoter Motifs
RAW_MOTIFS = {
    "TATA box":          "TATAWAWR",      # Classic TATA
    "Initiator (Inr)":   "YYANWYY",       # Initiator around TSS
    "DPE":               "RGWCGY",        # Downstream Promoter Element
    "BREu":              "SSRRCGCC",      # Upstream TFIIB Recognition Element
    "MTE":               "CSARCSS",       # Motif Ten Element
    "CAAT box":          "CCAAT",         # Proximal CAAT box
}

MOTIFS = {name: parse_iupac(seq) for name, seq in RAW_MOTIFS.items()}


def read_fasta(path):
    seqs = {}
    name = None
    if not os.path.exists(path):
        return seqs
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            else:
                seqs[name] += line.upper()
    return seqs


def shuffle_sequence_composition(seqs_dict):
    """Generates randomized background sequences matching length & base composition."""
    shuffled_dict = {}
    for name, seq in seqs_dict.items():
        chars = list(seq)
        random.shuffle(chars)
        shuffled_dict[name] = "".join(chars)
    return shuffled_dict


def scan_motifs_and_positions(seqs_dict):
    """Scans sequences for motifs on both strands, bins hits into the upstream-only windows,
    and logs individual hits."""
    results = {}
    all_hits_list = []

    for name, regex in MOTIFS.items():
        total_hits = plus_hits = minus_hits = 0
        pos_counts_plus = np.zeros(FULL_LEN)    # index 0 = -UPSTREAM_LEN ... FULL_LEN-1 = -1
        pos_counts_minus = np.zeros(FULL_LEN)
        win = {key: 0 for key, _, _, _ in WINDOWS}

        for seq_id, seq in seqs_dict.items():
            seq_len = len(seq)
            rev_seq = reverse_complement(seq)

            for strand, scan_seq in (("+", seq), ("-", rev_seq)):
                for m in regex.finditer(scan_seq):
                    match_len = m.end() - m.start()
                    if strand == "+":
                        start_idx = m.start()
                    else:
                        # map reverse-complement index back to forward-strand coordinates
                        start_idx = seq_len - m.start() - match_len
                    total_hits += 1
                    # TSS is the boundary right after the LAST base (last base = -1); anchoring
                    # on seq_len keeps promoters clipped shorter than 1000 bp correct.
                    rel_pos = start_idx - seq_len
                    dens_idx = rel_pos + UPSTREAM_LEN
                    if strand == "+":
                        plus_hits += 1
                        if 0 <= dens_idx < FULL_LEN:
                            pos_counts_plus[dens_idx] += 1
                    else:
                        minus_hits += 1
                        if 0 <= dens_idx < FULL_LEN:
                            pos_counts_minus[dens_idx] += 1

                    wkey = window_of(rel_pos)
                    if wkey:
                        win[wkey] += 1

                    all_hits_list.append({
                        "seq_id": seq_id,
                        "motif": name,
                        "strand": strand,
                        "start_pos_0based": start_idx,
                        "rel_tss_pos_bp": rel_pos,
                        "matched_sequence": m.group(0)
                    })

        results[name] = {
            'total_hits': total_hits,
            'plus_hits': plus_hits,
            'minus_hits': minus_hits,
            'win_distal': win["distal"],
            'win_proximal': win["proximal"],
            'win_core': win["core"],
            'pos_counts_plus': pos_counts_plus,
            'pos_counts_minus': pos_counts_minus,
            'pos_counts_total': pos_counts_plus + pos_counts_minus
        }

    return results, all_hits_list


SUMMARY_COLUMNS = [
    "motif", "real_total_hits", "rand_total_hits", "snr",
    "plus_strand_hits", "minus_strand_hits",
    f"distal_hits_{UPSTREAM_LEN}_to_151", "proximal_hits_150_to_31", "core_hits_30_to_1",
]


def build_summary_rows(real_res, rand_res):
    """One dict per motif with the SUMMARY_COLUMNS keys. Shared with 06b-organism."""
    rows = []
    for motif in MOTIFS:
        r, rnd = real_res[motif], rand_res[motif]
        snr = r['total_hits'] / max(1, rnd['total_hits'])
        rows.append({
            "motif": motif,
            "real_total_hits": r['total_hits'], "rand_total_hits": rnd['total_hits'],
            "snr": round(snr, 4),
            "plus_strand_hits": r['plus_hits'], "minus_strand_hits": r['minus_hits'],
            SUMMARY_COLUMNS[6]: r['win_distal'],
            SUMMARY_COLUMNS[7]: r['win_proximal'],
            SUMMARY_COLUMNS[8]: r['win_core'],
        })
    return rows


def plot_core_promoter_figure(real_res, num_promoters, label, out_path):
    """The 4-panel figure (strand split, upstream-only window split, TATA and Inr positional
    density). `label` prefixes the title (e.g. a genus name, or 'Pooled')."""
    motif_names = list(MOTIFS.keys())
    x = np.arange(len(motif_names))
    width = 0.25
    pre = f"{label} \u2014 " if label else ""

    fig, axes = plt.subplots(2, 2, figsize=(15, 11.5))
    fig.suptitle(f"{label + ' ' if label else ''}Core Promoter Element Analysis "
                 f"({num_promoters} Upstream-only Promoters)", fontsize=14, fontweight='bold')

    ax1 = axes[0, 0]
    ax1.bar(x - width / 2, [real_res[m]['plus_hits'] for m in motif_names], width,
            label='Plus (+) Strand', color='navy')
    ax1.bar(x + width / 2, [real_res[m]['minus_hits'] for m in motif_names], width,
            label='Minus (-) Strand', color='teal')
    ax1.set_ylabel('Motif Hits')
    ax1.set_title('A. Motif Hits by Strand (+ vs. -)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(motif_names, rotation=25, ha='right')
    ax1.legend(fontsize=9)
    ax1.grid(axis='y', linestyle='--', alpha=0.5)

    ax2 = axes[0, 1]
    colors = {"distal": "mediumpurple", "proximal": "darkorange", "core": "darkgreen"}
    for k, (key, _, _, wlabel) in enumerate(WINDOWS):
        ax2.bar(x + (k - 1) * width, [real_res[m]['win_' + key] for m in motif_names], width,
                label=wlabel, color=colors[key])
    ax2.set_ylabel('Motif Hits')
    ax2.set_title('B. Motif Distribution Across Upstream Sub-windows')
    ax2.set_xticks(x)
    ax2.set_xticklabels(motif_names, rotation=25, ha='right')
    ax2.legend(fontsize=8)
    ax2.grid(axis='y', linestyle='--', alpha=0.5)

    def plot_positional_profile(ax, motif_name, title):
        pos_total = real_res[motif_name]['pos_counts_total']
        pos_plus = real_res[motif_name]['pos_counts_plus']
        pos_minus = real_res[motif_name]['pos_counts_minus']
        offsets = np.arange(-UPSTREAM_LEN, DOWNSTREAM_LEN)
        if len(offsets) != len(pos_total):
            offsets = np.arange(len(pos_total)) - UPSTREAM_LEN
        w = 15
        k = np.ones(w) / w
        ax.plot(offsets, np.convolve(pos_total, k, mode='same'), 'k-', label='Total (+ and -)', linewidth=2.0)
        ax.plot(offsets, np.convolve(pos_plus, k, mode='same'), 'b--', label='Plus Strand (+)', linewidth=1.2)
        ax.plot(offsets, np.convolve(pos_minus, k, mode='same'), 'r--', label='Minus Strand (-)', linewidth=1.2)
        ax.axvspan(-30, 0, color='green', alpha=0.15, label='Core (-30 to -1)')
        ax.axvline(0, color='red', linestyle=':', linewidth=1.5, label='TSS (0)')
        ax.set_xlim(-UPSTREAM_LEN, 0)
        ax.set_xlabel('Position Relative to TSS (bp)')
        ax.set_ylabel('Motif Frequency (15-bp smoothed)')
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, linestyle=':', alpha=0.6)

    plot_positional_profile(axes[1, 0], 'TATA box', f'C. {pre}TATA Box Positional Density Profile')
    plot_positional_profile(axes[1, 1], 'Initiator (Inr)', f'D. {pre}Initiator (Inr) Positional Density Profile')

    fig.text(0.5, 0.005, UPSTREAM_ONLY_NOTE, ha='center', va='bottom', fontsize=8, style='italic', wrap=True)
    plt.tight_layout(rect=[0, 0.025, 1, 0.97])
    plt.savefig(out_path, dpi=300)
    plt.close()


def run_pooled_promoter_scan():
    real_seqs = read_fasta(INPUT_FASTA)

    if not real_seqs:
        print(f"Error: Could not find promoter sequences at '{INPUT_FASTA}'.")
        print("Please ensure the promoter-fetch stage has produced it.")
        return

    num_promoters = len(real_seqs)
    print(f"Loaded {num_promoters} upstream-only promoter sequences (pooled across all organisms).")

    rand_seqs = shuffle_sequence_composition(real_seqs)

    print("Scanning real promoters (+/- strands)...")
    real_res, real_hits = scan_motifs_and_positions(real_seqs)
    print("Scanning randomized background controls...")
    rand_res, _ = scan_motifs_and_positions(rand_seqs)

    print("\n" + "=" * 135)
    print(f"{'MOTIF':<18} | {'STRAND (+/-)':<14} | {'REAL':<6} | {'RAND':<6} | {'SNR':<6} | "
          f"{WINDOWS[0][3]:<22} | {WINDOWS[1][3]:<22} | {WINDOWS[2][3]:<16}")
    print("=" * 135)
    for motif in MOTIFS:
        r, rnd = real_res[motif], rand_res[motif]
        snr = r['total_hits'] / max(1, rnd['total_hits'])
        strand_str = f"+{r['plus_hits']} / -{r['minus_hits']}"
        print(f"{motif:<18} | {strand_str:<14} | {r['total_hits']:<6} | {rnd['total_hits']:<6} | "
              f"{snr:<6.2f} | {r['win_distal']:<22} | {r['win_proximal']:<22} | {r['win_core']:<16}")
        print("-" * 135)

    os.makedirs("outputs", exist_ok=True)

    with open(OUT_CSV_SUMMARY, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(build_summary_rows(real_res, rand_res))

    with open(OUT_CSV_HITS, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["seq_id", "motif", "strand", "start_pos_0based", "rel_tss_pos_bp", "matched_sequence"])
        for hit in real_hits:
            writer.writerow([hit["seq_id"], hit["motif"], hit["strand"],
                             hit["start_pos_0based"], hit["rel_tss_pos_bp"], hit["matched_sequence"]])

    print(f"\n[\u2713] Summary metrics exported to -> {OUT_CSV_SUMMARY}")
    print(f"[\u2713] Individual motif hits exported to -> {OUT_CSV_HITS}")

    plot_core_promoter_figure(real_res, num_promoters, "Pooled", OUT_PLOT)
    print(f"[\u2713] Analysis complete! Output figure saved to -> {OUT_PLOT}")


run_plasmodium_promoter_scan = run_pooled_promoter_scan   # backward-compatible alias


if __name__ == "__main__":
    run_pooled_promoter_scan()
