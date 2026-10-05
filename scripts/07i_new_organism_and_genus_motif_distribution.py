#!/usr/bin/env python3
"""
STAGE 16 (organism/genus) — Motif Distribution Block Diagram, built
separately for each organism AND each genus.

Follows the "separate image per organism" pattern (03b/04b, and 06b once
corrected) rather than a single comparative figure. A motif-block diagram
is fundamentally a per-sequence-set layout -- each row is one specific
promoter sequence -- so overlaying multiple organisms into one figure
wouldn't make sense here any more than it would for a tree or an MSA. Each
organism and each genus gets its own full diagram, identical in design to
07i's original single-figure output, just repeated per group.

REUSE: parse_full_meme_data, draw_motif_block, draw_scale_bar, and 07i's
own main() are all reused directly, unmodified. No multiprocessing here
(like 04b/06a/06b/07d, unlike 03b), so nothing needs duplicating. main() is
called once per group by redirecting its two module-level path constants
(MEME_XML_PATH, OUT_PLOT) before each call -- same technique 03b already
uses for stage 03's NUM_BOOTSTRAPS.

Falls back to a single "All_organisms" legacy run if neither a per-organism
nor a per-genus meme.xml split exists, but a bare outputs/meme_output/meme.xml
does -- same fallback convention as the GOMO stage (07d).

KNOWN MINOR LIMITATION: the in-image plot title doesn't say which
organism/genus it is -- 07i's main() hardcodes that title text with no
label parameter, and this script leaves 07i completely unmodified,
matching the "never touch the original" precedent from every other stage
in this project. Rely on the filename (<Label>_motif_distribution.png)
for that instead.

Requirements:
    pip install numpy matplotlib seaborn
"""

import os
import re
import sys
import glob
import importlib.util

# ============================== CONFIGURATION ==============================
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

STAGE07I_FILENAME = None   # set explicitly to skip auto-detection, e.g.
                            # "07i_motif_distribution_block_diagram.py"

MEME_OUT_DIR = "outputs/meme_output"                 # species-level (07c)
MEME_OUT_DIR_GENUS = "outputs/meme_output_genus"      # genus-level (07c)
LEGACY_MEME_XML = "outputs/meme_output/meme.xml"      # pre-split fallback

PLOT_DIR_ORGANISM = "outputs/motif_distribution_by_organism"
PLOT_DIR_GENUS = "outputs/motif_distribution_by_genus"
# ===========================================================================


def _find_stage07i_module():
    if STAGE07I_FILENAME:
        path = os.path.join(_THIS_DIR, STAGE07I_FILENAME)
        if not os.path.exists(path):
            sys.exit(f"STAGE07I_FILENAME is set to {STAGE07I_FILENAME!r} but that "
                      f"file doesn't exist in {_THIS_DIR}.")
        return path

    candidates = sorted(
        p for p in glob.glob(os.path.join(_THIS_DIR, "07i_*.py"))
        if os.path.basename(p) != os.path.basename(__file__)
    )
    if not candidates:
        sys.exit(f"Could not auto-detect the stage-07i motif-distribution script "
                  f"(07i_*.py) in {_THIS_DIR}. Set STAGE07I_FILENAME above.")
    if len(candidates) > 1:
        print("Multiple 07i_*.py scripts found:")
        for p in candidates:
            print(f"    {os.path.basename(p)}")
        print(f"Using {os.path.basename(candidates[0])}. Set STAGE07I_FILENAME above to override.\n")
    return candidates[0]


_STAGE07I_PATH = _find_stage07i_module()
_spec = importlib.util.spec_from_file_location("motif_dist_base", _STAGE07I_PATH)
mdb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mdb)  # only defines functions/constants -- its __main__
                                # guard keeps its own legacy single-file run from
                                # happening here as an import side effect.


def discover_meme_xml_groups(meme_out_dir):
    groups = {}
    for xml_path in sorted(glob.glob(os.path.join(meme_out_dir, "*", "meme.xml"))):
        label = os.path.basename(os.path.dirname(xml_path))
        groups[label] = xml_path
    return groups


def sanitize_filename(name):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip())
    return safe.strip("_") or "Unknown"


def run_for_group(label, meme_xml_path, out_path, tag_word):
    print(f"\n{'=' * 60}\n[{tag_word.capitalize()}] {label}\n{'=' * 60}")
    # 07i's own main() assumes OUT_PLOT sits directly under a flat "outputs/"
    # (it only ever calls os.makedirs("outputs", ...)) -- since this script
    # points it at a nested per-group subfolder instead, that subfolder is
    # created here first, or plt.savefig() would fail with FileNotFoundError.
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    mdb.MEME_XML_PATH = meme_xml_path
    mdb.OUT_PLOT = out_path
    try:
        mdb.main()
        return True
    except SystemExit as e:
        print(f"    [!] Skipped {label} (missing/empty data): {e}")
        return False
    except Exception as e:
        print(f"    [!] Skipped {label} (unexpected error): {e}")
        return False


def main():
    species_groups = discover_meme_xml_groups(MEME_OUT_DIR)
    genus_groups = discover_meme_xml_groups(MEME_OUT_DIR_GENUS)

    if not species_groups and not genus_groups:
        if os.path.exists(LEGACY_MEME_XML):
            print(f"[!] No per-organism/per-genus meme.xml found under {MEME_OUT_DIR}*/ -- "
                  f"falling back to the single legacy file at {LEGACY_MEME_XML}.")
            species_groups = {"All_organisms": LEGACY_MEME_XML}
        else:
            sys.exit(f"[!] No meme.xml files found under '{MEME_OUT_DIR}/<Organism>/', "
                      f"'{MEME_OUT_DIR_GENUS}/<Genus>/', or '{LEGACY_MEME_XML}'. "
                      f"Run 07c_Motif_discovery_for_all.py first.")

    print(f"[+] Found {len(species_groups)} organism-level and {len(genus_groups)} "
          f"genus-level meme.xml file(s).")

    print(f"\n{'#' * 60}\n# SPECIES-LEVEL DIAGRAMS\n{'#' * 60}")
    n_done = 0
    for label in sorted(species_groups.keys()):
        out_path = os.path.join(PLOT_DIR_ORGANISM, f"{sanitize_filename(label)}_motif_distribution.png")
        n_done += int(run_for_group(label, species_groups[label], out_path, "organism"))
    print(f"\n[\u2713] {n_done}/{len(species_groups)} organism diagram(s) -> {PLOT_DIR_ORGANISM}/")

    if genus_groups:
        print(f"\n{'#' * 60}\n# GENUS-LEVEL DIAGRAMS\n{'#' * 60}")
        n_done_genus = 0
        for label in sorted(genus_groups.keys()):
            out_path = os.path.join(PLOT_DIR_GENUS, f"{sanitize_filename(label)}_motif_distribution.png")
            n_done_genus += int(run_for_group(label, genus_groups[label], out_path, "genus"))
        print(f"\n[\u2713] {n_done_genus}/{len(genus_groups)} genus diagram(s) -> {PLOT_DIR_GENUS}/")
    else:
        print("\n[i] No genus-level meme.xml files found -- skipped.")


if __name__ == "__main__":
    main()
