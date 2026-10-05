#!/usr/bin/env python3
"""
STAGE 07i (v3) — Motif Distribution Block Diagrams, per ORGANISM and per GENUS.

WHAT CHANGED vs. v2
---------------------------------------------------------------------------
1. Layout. 07c no longer writes a single outputs/meme_output/meme.xml. It
   writes one meme.xml per organism  (outputs/meme_output/<Organism>/) and
   one per genus (outputs/meme_output_genus/<Genus>/). This script now
   auto-discovers every one of them (same discover_meme_xml_groups() approach
   as 07e / 07f) and draws one diagram per group:

     outputs/motif_distribution_by_organism/<Organism>_motif_distribution.png
     outputs/motif_distribution_by_genus/<Genus>_motif_distribution.png

   If neither layout exists, it falls back to the legacy single
   outputs/meme_output/meme.xml (drawn as group "All_organisms").

2. TSS-relative coordinates. MEME's `position` is a 1-based offset from the
   START of each (upstream-only) promoter record. The diagram now plots every
   site TSS-relative, from the sequence's own length in meme.xml:

       site_start = position - seq_len - 1        (last base of the record = -1)

   so the x-axis runs from about -1000 (left) to the TSS at 0 (right), exactly
   like 06a/06b/06c/07f, and promoters clipped shorter than 1000 bp are
   correctly right-aligned on the TSS instead of being drawn as if they
   started at the same place.

3. Strand arrows are meaningful now that 07c runs MEME with -revcomp: `plus`
   sites point right, `minus` sites point left (MEME writes strand="plus" /
   "minus" in meme.xml; "+" / "-" are accepted too). Sites are assumed to be
   reported by their leftmost forward-strand coordinate for both strands.

4. The custom ruler is the ONLY x-axis (the default matplotlib axis is
   disabled), as in v2.

Requirements:  pip install numpy matplotlib seaborn
"""

import os
import re
import sys
import glob
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

# ============================== CONFIGURATION ==============================
MEME_OUT_DIR = "outputs/meme_output"                 # per-organism (07c)
MEME_OUT_DIR_GENUS = "outputs/meme_output_genus"      # per-genus (07c)
LEGACY_MEME_XML = "outputs/meme_output/meme.xml"      # pre-split fallback

OUT_DIR_ORG = "outputs/motif_distribution_by_organism"
OUT_DIR_GENUS = "outputs/motif_distribution_by_genus"

DEFAULT_PROMOTER_LENGTH = 1000    # used only if a <sequence> has no length attribute
RULER_TICK_STEP = 100
MAX_FIG_HEIGHT_IN = 60            # cap for very large groups
# ===========================================================================


def sanitize_filename(name):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip())
    return safe.strip("_") or "Unknown"


def discover_meme_xml_groups(meme_out_dir):
    """Finds every <meme_out_dir>/<Label>/meme.xml and returns {label: path}."""
    groups = {}
    for xml_path in sorted(glob.glob(os.path.join(meme_out_dir, "*", "meme.xml"))):
        groups[os.path.basename(os.path.dirname(xml_path))] = xml_path
    return groups


def parse_full_meme_data(xml_path):
    """Returns (seq_info, motifs_meta, occurrences) or None if unusable.
    Every occurrence carries TSS-relative start/end (see module docstring)."""
    try:
        root = ET.parse(xml_path).getroot()
    except Exception as e:
        print(f"    [!] Failed to parse {xml_path}: {e}")
        return None

    seq_info = {}
    for seq in root.findall(".//training_set/sequence"):
        try:
            length = int(seq.get("length"))
        except (TypeError, ValueError):
            length = DEFAULT_PROMOTER_LENGTH
        seq_info[seq.get("id")] = {"name": seq.get("name", seq.get("id")), "length": length}
    if not seq_info:
        print(f"    [!] No <sequence> entries found in {xml_path}.")
        return None

    motifs_meta, occurrences = {}, []
    for motif in root.findall(".//motifs/motif"):
        m_id = motif.get("id")
        width = int(motif.get("width", 0))
        motifs_meta[m_id] = {
            "name": motif.get("name", m_id),
            "width": width,
            "consensus": motif.get("alt") or motif.get("name", m_id),
            "e_value": motif.get("e_value", ""),
        }
        for site in motif.findall(".//contributing_sites/contributing_site"):
            sid = site.get("sequence_id")
            if sid not in seq_info:
                continue
            seq_len = seq_info[sid]["length"]
            raw_pos = int(site.get("position"))
            start_rel = raw_pos - seq_len - 1          # last base of record -> -1
            strand_raw = (site.get("strand") or "plus").lower()
            strand = "-" if strand_raw in ("minus", "-") else "+"
            occurrences.append({
                "seq_id": sid, "motif_id": m_id, "width": width,
                "raw_position": raw_pos,
                "start_rel": start_rel, "end_rel": start_rel + width - 1,
                "strand": strand,
            })

    if not motifs_meta:
        print(f"    [!] No <motif> entries found in {xml_path}.")
        return None
    return seq_info, motifs_meta, occurrences


def draw_motif_block(ax, y, start, width, strand, color, edge_frac=0.2):
    """Block spanning [start, start+width] on the x axis (continuous coordinates; base -1
    occupies [-1, 0]). The arrow head sits at the right end for '+' and the left end for '-'."""
    height = 0.6
    head = max(width * edge_frac, 0.5)
    body = width - head
    if strand == "+":
        rect_x = start
        tri_x = [start + body, start + body, start + width]
    else:
        rect_x = start + head
        tri_x = [start + head, start + head, start]
    ax.add_patch(patches.Rectangle((rect_x, y - height / 2), body, height,
                                   facecolor=color, edgecolor="black", linewidth=0.4, zorder=3))
    ax.add_patch(patches.Polygon(list(zip(tri_x, [y - height / 2, y + height / 2, y])), closed=True,
                                 facecolor=color, edgecolor="black", linewidth=0.4, zorder=3))


def draw_ruler(ax, y_bottom, x_min, x_max, tick_step):
    """The ONLY x-axis in the figure: TSS-relative bp, ticks every tick_step, 0 = TSS."""
    ax.plot([x_min, x_max], [y_bottom, y_bottom], color="black", linewidth=1.5, zorder=4)
    first_tick = -((-x_min) // tick_step) * tick_step   # first multiple of tick_step >= x_min
    for tick in range(int(first_tick), int(x_max) + 1, tick_step):
        ax.plot([tick, tick], [y_bottom - 0.15, y_bottom + 0.15], color="black", linewidth=1.2, zorder=4)
        ax.text(tick, y_bottom - 0.45, "TSS" if tick == 0 else str(tick),
                ha="center", va="top", fontsize=7, zorder=4)
    ax.text((x_min + x_max) / 2, y_bottom - 1.15,
            "Position relative to TSS (bp; -1 = base immediately upstream)",
            ha="center", va="top", fontsize=9, fontweight="bold", zorder=4)


def plot_group(label, xml_path, out_dir, tag_word):
    print(f"\n{'-' * 60}\n{tag_word.capitalize()}: {label}\n{'-' * 60}")
    parsed = parse_full_meme_data(xml_path)
    if parsed is None:
        return False
    seq_info, motifs_meta, occurrences = parsed

    seq_ids = list(seq_info.keys())
    n = len(seq_ids)
    n_minus = sum(1 for o in occurrences if o["strand"] == "-")
    print(f"    [*] {n} sequences, {len(motifs_meta)} motifs, {len(occurrences)} sites "
          f"({len(occurrences) - n_minus} plus / {n_minus} minus).")
    if not occurrences:
        print("    [!] MEME reported no contributing sites -- skipping.")
        return False

    bad = [o for o in occurrences
           if o["end_rel"] > -1 or o["start_rel"] < -seq_info[o["seq_id"]]["length"]]
    if bad:
        print(f"    [!] {len(bad)} site(s) fall outside [-seq_len, -1] after conversion -- "
              f"check that meme.xml was built from the upstream-only promoter FASTA.")

    motif_ids = sorted(motifs_meta.keys())
    palette = sns.color_palette("Set2", n_colors=max(3, len(motif_ids)))
    color_map = dict(zip(motif_ids, palette))

    occ_by_seq = {}
    for o in occurrences:
        occ_by_seq.setdefault(o["seq_id"], []).append(o)

    max_len = max(info["length"] for info in seq_info.values())
    fig_h = min(MAX_FIG_HEIGHT_IN, max(6, n * 0.22 + 2.5))
    fig, ax = plt.subplots(figsize=(14, fig_h), dpi=200)

    for row, sid in enumerate(seq_ids):
        y = n - row
        L = seq_info[sid]["length"]
        ax.plot([-L, 0], [y, y], color="#cccccc", linewidth=3, zorder=1, solid_capstyle="butt")
        for o in occ_by_seq.get(sid, []):
            draw_motif_block(ax, y, o["start_rel"], o["width"], o["strand"], color_map[o["motif_id"]])

    ax.axvline(0, color="red", linestyle=":", linewidth=1.2, zorder=2)
    ax.set_yticks(range(1, n + 1))
    ax.set_yticklabels([seq_info[s]["name"] for s in reversed(seq_ids)],
                       fontsize=7 if n <= 60 else 5)
    ax.set_xlim(-max_len - 30, 30)
    ax.set_ylim(-2.0, n + 1.5)
    ax.set_xticks([])
    for side in ("bottom", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title(f"Motif Distribution \u2014 {tag_word.capitalize()}: {label} "
                 f"({n} sequences, {len(motif_ids)} motifs; arrows = strand)",
                 fontsize=12, fontweight="bold", pad=12)

    draw_ruler(ax, y_bottom=-1.0, x_min=-max_len, x_max=0, tick_step=RULER_TICK_STEP)

    handles = [patches.Patch(facecolor=color_map[m], edgecolor="black",
                             label=f"{m} ({motifs_meta[m]['consensus']})") for m in motif_ids]
    ax.legend(handles=handles, title="Motif", loc="upper left", fontsize=7,
              framealpha=0.9, bbox_to_anchor=(1.01, 1.0))

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sanitize_filename(label)}_motif_distribution.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    [\u2713] Diagram saved -> {out_path}")
    return True


def main():
    os.makedirs("outputs", exist_ok=True)
    org_groups = discover_meme_xml_groups(MEME_OUT_DIR)
    genus_groups = discover_meme_xml_groups(MEME_OUT_DIR_GENUS)

    if not org_groups and not genus_groups:
        if os.path.exists(LEGACY_MEME_XML):
            print(f"[!] No per-organism/per-genus meme.xml found -- falling back to {LEGACY_MEME_XML}.")
            org_groups = {"All_organisms": LEGACY_MEME_XML}
        else:
            sys.exit(f"[!] No meme.xml found under '{MEME_OUT_DIR}/<Organism>/', "
                     f"'{MEME_OUT_DIR_GENUS}/<Genus>/' or '{LEGACY_MEME_XML}'. Run 07c first.")

    print(f"[+] Found {len(org_groups)} organism-level and {len(genus_groups)} genus-level meme.xml file(s).")

    print(f"\n{'#' * 60}\n# ORGANISM-LEVEL DIAGRAMS\n{'#' * 60}")
    ok_o = sum(plot_group(l, p, OUT_DIR_ORG, "organism") for l, p in sorted(org_groups.items()))
    print(f"\n{'#' * 60}\n# GENUS-LEVEL DIAGRAMS\n{'#' * 60}")
    ok_g = sum(plot_group(l, p, OUT_DIR_GENUS, "genus") for l, p in sorted(genus_groups.items()))

    print(f"\n[\u2713] {ok_o}/{len(org_groups)} organism-level and {ok_g}/{len(genus_groups)} "
          f"genus-level diagram(s) written.")
    print(f"    -> {OUT_DIR_ORG}/<Organism>_motif_distribution.png")
    print(f"    -> {OUT_DIR_GENUS}/<Genus>_motif_distribution.png")


if __name__ == "__main__":
    main()
