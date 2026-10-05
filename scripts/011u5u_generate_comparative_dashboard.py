#!/usr/bin/env python3
"""
STAGE 011c — Comparative Dashboard (genus-wise by default, species-wise
available as an on-the-fly download).

WHAT CHANGED FROM THE PREVIOUS TWO-FILE VERSION
---------------------------------------------------------------------------
This used to write TWO separate HTML files (genus + species), each linking
to the other for download. It now writes exactly ONE file:
    promoter_comparative_dashboard.html
which shows the genus-wise comparison by default. Its top-right toolbar has
a single button, "Download species-wise dashboard" -- clicking it builds
and downloads a complete, standalone species-wise HTML file entirely in
the browser (via a JS Blob + object URL), with NO second file ever
written to disk by this script.

HOW "ON THE FLY, NO SECOND FILE" ACTUALLY WORKS (read this before editing)
---------------------------------------------------------------------------
This is a static HTML dashboard with no live backend -- once the browser
has the page open, this Python script isn't running anymore, so nothing
can be freshly re-scanned from outputs/ at click time. "On the fly" here
means: at BUILD time, this script still renders the complete species-wise
document (same panels/images/lightbox as the old species file) exactly as
before, but instead of writing it to its own file, embeds that finished
HTML as a single JS string constant (SPECIES_DASHBOARD_HTML) inside the
ONE file that does get written. The download button just hands that
already-built string to the browser as a Blob and triggers a save -- so
disk-wise there is only ever one file, but the generation cost of the
species view is paid once, up front, same as before.

CRITICAL SUBTLETY: escaping "</script" inside the embedded string
---------------------------------------------------------------------------
The embedded species document contains its own <script>...</script> block
(the lightbox JS). If that literal "</script>" text were embedded verbatim
inside the OUTER page's <script> tag, the browser's HTML tokenizer -- which
looks for that exact substring textually, with no awareness of JS string
literals -- would think the OUTER script tag ends right there, silently
truncating everything after it. json.dumps() alone does NOT protect
against this (JSON escaping has no rule for "</script"). So after
json.dumps()-encoding the species HTML, every occurrence of the raw
substring "</script" is replaced with "<\\/script" (an escaped forward
slash -- JS treats a backslash-escaped slash identically to a plain slash, so the string's *meaning* is
unchanged, but the literal textual match the HTML parser looks for is
broken). This is applied generically to the whole encoded string, so it
stays correct even if the embedded document's own markup changes later.

Everything else -- panel definitions, candidate path patterns, genus/species
auto-detection from family_labels.csv, the CSS theme, the lightbox -- is
unchanged from the previous version. See the per-panel comments below for
which filenames are confirmed vs. educated-guess globs.

Requirements: none beyond the standard library.
"""

import os
import re
import csv
import glob
import json
import base64
import xml.etree.ElementTree as ET
from pathlib import Path

OUTPUT_DIR = "./outputs"
LABELS_PATH = os.path.join(OUTPUT_DIR, "family_labels.csv")
PROMOTER_FASTA_PATH = os.path.join(OUTPUT_DIR, "promoter_sequences.fasta")

MEME_OUTPUT_DIR_ORGANISM = os.path.join(OUTPUT_DIR, "meme_output")
MEME_OUTPUT_DIR_GENUS = os.path.join(OUTPUT_DIR, "meme_output_genus")
TOMTOM_OUTPUT_DIR_ORGANISM = os.path.join(OUTPUT_DIR, "tomtom_output")
TOMTOM_OUTPUT_DIR_GENUS = os.path.join(OUTPUT_DIR, "tomtom_output_genus")
CPG_DETAILED_DIR_ORGANISM = os.path.join(OUTPUT_DIR, "CpG_Island_By_Organism")

# A motif's best TOMTOM hit is only treated as "matched to a known TF" (and
# labeled/colored by that TF's name in the Architecture Map) if its q-value
# is at or below this. Above it, the motif is shown as "Motif N (unmatched)"
# instead -- a real hit exists in the TSV either way, this only controls
# the significance bar for LABELING it as a confident match.
TOMTOM_QVALUE_CUTOFF = 0.05

# The ONE file this script writes to disk.
DASHBOARD_FILE = "promoter_comparative_dashboard.html"
# Filename suggested to the browser when the download button fires --
# this is never written by this script itself, only ever produced
# client-side from the embedded SPECIES_DASHBOARD_HTML string.
SPECIES_DOWNLOAD_FILENAME = "promoter_comparative_dashboard_species.html"

# Set to an explicit list to override auto-detection. None = auto-detect
# every genus (GENERA_OVERRIDE) / every full organism name
# (SPECIES_OVERRIDE) present in outputs/family_labels.csv's "organism"
# column.
GENERA_OVERRIDE = None
SPECIES_OVERRIDE = None

# If True, a genus-dashboard panel that only has organism-level (not
# genus-level) data shows EVERY species under that genus side by side,
# instead of just the first one found. Off by default to keep panels a
# manageable size. Has no effect on the species-wise dashboard.
SHOW_ALL_SPECIES_PER_GENUS = False

# ============================================================================
# Panel definitions. Each panel: a label, and an ordered list of glob
# PATTERNS (relative to OUTPUT_DIR, "{group}" is substituted per genus/
# species). The FIRST pattern that matches anything for a given group is
# used.
# ============================================================================

GENUS_SECTIONS = [
    {
        "title": "1. Protein Physicochemical Properties",
        "panels": [
            {
                "label": "Protein Property Distributions (genus-wise)",
                # CONFIRMED exact filename/folder -- directly from your own
                # `ls protein_properties_by_genus/Toxoplasma/` output.
                "patterns": [
                    "protein_properties_by_genus/{group}/protein_properties_distributions.png",
                ],
            },
        ],
    },
    {
        "title": "2. Evolutionary & Conservation",
        "panels": [
            {
                "label": "MSA Conservation Profile (genus-wise)",
                # CONFIRMED exact folder/filename from 04b's own log output:
                # "phylo_by_organism/<Genus>/<Genus>_msa_conservation_zoom.png"
                # (folder is misleadingly named "*_by_organism" but 04b's log
                # confirms it holds one alignment PER GENUS, not per species).
                # This is listed FIRST and is intentionally strict/exact --
                # the old loose fallback below ("**/{group}*msa*conservation*.png")
                # also matches any per-species file whose name starts with this
                # genus (e.g. "Toxoplasma_gondii_FOU_msa_conservation_zoom.png"
                # literally starts with "Toxoplasma"), and since only the first
                # matching pattern's (alphabetically) sorted first result is
                # used, that loose pattern could silently return a
                # species-level image instead of the genus one. Keeping the
                # exact pattern first means the fallback is only ever reached
                # if this confirmed path doesn't exist.
                "patterns": [
                    "phylo_by_organism/{group}/{group}_msa_conservation_zoom.png",
                    "phylo_by_organism/{group}/{group}_conservation_profile.png",
                    "**/*msa*by_genus*/{group}*conservation*.png",
                    "**/{group}*msa*conservation*.png",
                    "**/{group}*conservation*profile*.png",
                ],
            },
        ],
    },
    {
        "title": "3. Architecture & Synergism",
        "panels": [
            {
                "label": "Core Promoter Element Analysis",
                # New panel -- Core_Promoter_Results_By_Organism/ existed in
                # your outputs but had no panel at all before. Same
                # genus-pooled-folder situation as Synergism below: your
                # listing shows only a bare "Toxoplasma" subfolder (not a
                # per-species one) under this "_By_Organism"-named path, so
                # the genus-wise glob (which just needs {group}* to match
                # that folder) already works correctly here.
                "patterns": [
                    "Core_Promoter_Results_By_Genus/{group}*/**/*.png",
                    "Core_Promoter_Results_By_Organism/{group}*/**/*.png",
                    "**/{group}*core*promoter*result*.png",
                ],
            },
            {
                "label": "Core Promoter Element Synergism",
                # Confirmed from your listing: Core_Promoter_Synergism_By_Organism
                # contains a bare "Toxoplasma" folder (genus-pooled data,
                # despite the "_By_Organism" name) -- "{group}*" with
                # group="Toxoplasma" matches that folder exactly, so this
                # already works for the genus dashboard as-is.
                "patterns": [
                    "Core_Promoter_Synergism_By_Genus/{group}*/**/*.png",
                    "Core_Promoter_Synergism_By_Organism/{group}*/**/*.png",
                    "**/{group}*synergism*.png",
                ],
            },
            {
                "label": "CpG Island Distribution (genus-wise)",
                # Confirmed exact filename + folder from
                # 08_new_organism_wise_cpg_island_finder.py.
                "patterns": [
                    "CpG_Island_By_Genus/{group}/{group}_cpg_island_distribution.png",
                    "**/{group}_cpg_island_distribution.png",
                ],
            },
            {
                "label": "DNA Structural Properties (genus-wise)",
                "patterns": [
                    "DNA_Structural_Properties_By_Genus/{group}/**/*.png",
                    "**/{group}*structural*propert*.png",
                ],
            },
        ],
    },
    {
        "title": "4. Motif Discovery & Functional",
        "panels": [
            {
                "label": "Enriched Motif Logos (genus-wise)",
                # CORRECTNESS FIX: your listing shows
                # motif_logos_by_genus/Toxoplasma_logos.png is a FLAT FILE,
                # not a subfolder -- the old "{group}/**/*.png" pattern
                # (expects a subfolder) never matched it. This was wrong
                # despite being labeled "correct" before; fixed here using
                # your confirmed real filename, with the old subfolder
                # guess kept as a harmless secondary fallback in case a
                # future run nests differently.
                "patterns": [
                    "motif_logos_by_genus/{group}_logos.png",
                    "motif_logos_by_genus/{group}/**/*.png",
                    "**/{group}*motif*logo*.png",
                ],
            },
            {
                "label": "Motif Distribution Across Sequences (genus-wise)",
                # CORRECTNESS FIX: confirmed flat file
                # motif_distribution_by_genus/Toxoplasma_motif_distribution.png
                # in your listing -- same issue and same fix as Motif Logos above.
                "patterns": [
                    "motif_distribution_by_genus/{group}_motif_distribution.png",
                    "motif_distribution_by_genus/{group}/**/*.png",
                    "**/{group}*motif*distribution*.png",
                ],
            },
            {
                "label": "Motif Conservation Analysis (genus-wise)",
                # Not directly confirmed at the genus level in your ls
                # output (you didn't list inside motif_conservation_by_genus/),
                # but motif_conservation_by_organism/ is CONFIRMED flat
                # (see species section below) and both are written by the
                # same script with the same naming convention, so the flat
                # form is listed first with the old subfolder guess kept as
                # a fallback.
                "patterns": [
                    "motif_conservation_by_genus/{group}_motif_conservation_analysis.png",
                    "motif_conservation_by_genus/{group}/**/*.png",
                    "**/{group}*motif*conservation*.png",
                ],
            },
        ],
    },
]

# Species-wise mirrors the genus one, minus the phylogenetic tree panel
# (see module docstring), pointed at "_By_Organism" folders directly.
SPECIES_SECTIONS = [
    {
        "title": "1. Protein Physicochemical Properties",
        "panels": [
            {
                "label": "Protein Property Distributions (species-wise)",
                # CONFIRMED exact filename/folder -- directly from your own
                # `ls protein_properties_by_organism/Plasmodium_berghei_ANKA/`
                # output.
                "patterns": [
                    "protein_properties_by_organism/{group}/protein_properties_distributions.png",
                ],
            },
        ],
    },
    {
        "title": "2. Evolutionary & Conservation",
        "panels": [
            {
                "label": "MSA Conservation Profile (species-wise)",
                "patterns": [
                    "**/*msa*by_organism*/{group}*conservation*.png",
                    "**/{group}*msa*conservation*.png",
                    "**/{group}*conservation*profile*.png",
                ],
            },
        ],
    },
    {
        "title": "3. Architecture & Synergism",
        "panels": [
            {
                "label": "Core Promoter Element Analysis",
                # Your listing shows Core_Promoter_Results_By_Organism/
                # contains only a bare "Toxoplasma" folder -- despite its
                # "_By_Organism" name, this data is currently genus-pooled,
                # not actually per-species. The exact-species pattern is
                # tried first (future-proofs this if that stage is ever
                # fixed to be truly per-organism), then falls back to the
                # genus-pooled folder via {genus} so the panel shows
                # SOMETHING rather than "Not found" for every species --
                # clearly the same genus-pooled image repeated across every
                # species in that genus, not a real per-species result.
                "patterns": [
                    "Core_Promoter_Results_By_Organism/{group}/**/*.png",
                    "Core_Promoter_Results_By_Genus/{genus}/**/*.png",
                    "Core_Promoter_Results_By_Organism/{genus}/**/*.png",
                    "**/{group}*core*promoter*result*.png",
                ],
            },
            {
                "label": "Core Promoter Element Synergism",
                # Same genus-pooled-data situation as above, confirmed by
                # your listing (Core_Promoter_Synergism_By_Organism/ also
                # contains only a bare "Toxoplasma" folder, not
                # "Toxoplasma_gondii_ARI" etc.) -- same exact-species-first,
                # genus-pooled-fallback fix.
                "patterns": [
                    "Core_Promoter_Synergism_By_Organism/{group}/**/*.png",
                    "Core_Promoter_Synergism_By_Genus/{genus}/**/*.png",
                    "Core_Promoter_Synergism_By_Organism/{genus}/**/*.png",
                    "**/{group}*synergism*.png",
                ],
            },
            {
                "label": "CpG Island Distribution (species-wise)",
                # Confirmed exact filename pattern AND confirmed genuinely
                # per-species folder ("CpG_Island_By_Organism/Toxoplasma_gondii_ARI/"
                # etc.) directly from your listing.
                "patterns": [
                    "CpG_Island_By_Organism/{group}/{group}_cpg_island_distribution.png",
                    "**/{group}_cpg_island_distribution.png",
                ],
            },
            {
                "label": "DNA Structural Properties (species-wise)",
                # Confirmed genuinely per-species folder
                # ("DNA_Structural_Properties_By_Organism/Toxoplasma_gondii_ARI/")
                # from your listing; exact filename inside not confirmed,
                # so this stays a recursive glob.
                "patterns": [
                    "DNA_Structural_Properties_By_Organism/{group}/**/*.png",
                    "**/{group}*structural*propert*.png",
                ],
            },
        ],
    },
    {
        "title": "4. Motif Discovery & Functional",
        "panels": [
            {
                "label": "Enriched Motif Logos (species-wise)",
                # CORRECTNESS FIX: confirmed flat file
                # motif_logos_by_organism/Toxoplasma_gondii_ARI_logos.png
                # in your listing -- the old "{group}/**/*.png" pattern
                # (expects a subfolder) never matched this.
                "patterns": [
                    "motif_logos_by_organism/{group}_logos.png",
                    "motif_logos_by_organism/{group}/**/*.png",
                    "**/{group}*motif*logo*.png",
                ],
            },
            {
                "label": "Motif Distribution Across Sequences (species-wise)",
                # CORRECTNESS FIX: confirmed flat file
                # motif_distribution_by_organism/Toxoplasma_gondii_ARI_motif_distribution.png
                # in your listing -- same issue, same fix.
                "patterns": [
                    "motif_distribution_by_organism/{group}_motif_distribution.png",
                    "motif_distribution_by_organism/{group}/**/*.png",
                    "**/{group}*motif*distribution*.png",
                ],
            },
            {
                "label": "Motif Conservation Analysis (species-wise)",
                # CORRECTNESS FIX: confirmed flat file
                # motif_conservation_by_organism/Toxoplasma_gondii_ARI_motif_conservation_analysis.png
                # in your listing -- same issue, same fix.
                "patterns": [
                    "motif_conservation_by_organism/{group}_motif_conservation_analysis.png",
                    "motif_conservation_by_organism/{group}/**/*.png",
                    "**/{group}*motif*conservation*.png",
                ],
            },
        ],
    },
]

# ===========================================================================


def sanitize_organism_name(name):
    """Same alnum-only sanitizer used throughout the pipeline (01, 06b/06c,
    08-organism, 03c, 23_comparative_matrix.py) so this matches the actual
    folder/file names those scripts wrote."""
    if not name or str(name).strip() == "" or str(name).lower() == "nan":
        return "Unknown_organism"
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_")


def detect_genera(labels_path):
    """Auto-detects genera from family_labels.csv's 'organism' column
    (first word of each entry)."""
    if not os.path.exists(labels_path):
        print(f"[!] {labels_path} not found -- can't auto-detect genera.")
        return []
    genera = set()
    with open(labels_path) as f:
        for row in csv.DictReader(f):
            organism = (row.get("organism") or "").strip()
            if organism:
                genera.add(organism.split()[0])
    return sorted(genera)


def detect_species(labels_path):
    """Auto-detects every full organism name from family_labels.csv's
    'organism' column, sanitized the same way every organism-wise script's
    output folders already are."""
    if not os.path.exists(labels_path):
        print(f"[!] {labels_path} not found -- can't auto-detect species.")
        return []
    species = set()
    with open(labels_path) as f:
        for row in csv.DictReader(f):
            organism = (row.get("organism") or "").strip()
            if organism:
                species.add(sanitize_organism_name(organism))
    return sorted(species)


def find_images_for_group(group, patterns, search_dir, show_all):
    """Tries each pattern (in order) with {group} (and {genus}, the first
    underscore-separated token of group -- meaningful for species IDs like
    "Toxoplasma_gondii_ARI" -> "Toxoplasma"; for a genus group like
    "Toxoplasma" it's simply the same value) substituted, using recursive
    glob relative to search_dir. Returns a list of matched paths -- one
    file, unless show_all is True and the matching pattern's glob finds
    several."""
    genus = str(group).split("_")[0] if group else group
    for pattern in patterns:
        resolved = pattern.format(group=group, genus=genus)
        matches = sorted(glob.glob(os.path.join(search_dir, resolved), recursive=True))
        matches = [m for m in matches if os.path.isfile(m)]
        if matches:
            return matches if show_all else matches[:1]
    return []


def image_to_base64(img_path):
    try:
        ext = Path(img_path).suffix.lower().replace(".", "")
        if ext == "svg":
            ext = "svg+xml"
        with open(img_path, "rb") as img_file:
            encoded = base64.b64encode(img_file.read()).decode("utf-8")
            return f"data:image/{ext};base64,{encoded}"
    except Exception as e:
        print(f"Error reading image {img_path}: {e}")
        return None


def render_group_card(group_label, matches):
    if not matches:
        return f"""
            <div class="group-card missing-card">
                <div class="group-card-label">{group_label}</div>
                <div class="missing-msg">Not found</div>
            </div>"""

    cards = []
    for img_path in matches:
        b64_str = image_to_base64(img_path)
        if not b64_str:
            continue
        filename = os.path.basename(img_path)
        js_safe_filename = filename.replace("\\", "\\\\").replace("'", "\\'")
        cards.append(f"""
            <div class="group-card">
                <div class="group-card-label">{group_label}</div>
                <img src="{b64_str}" alt="{group_label}"
                     onclick="openLightbox('{b64_str}', '{js_safe_filename}')"/>
                <div class="card-caption">{filename}</div>
            </div>""")
    return "".join(cards) if cards else render_group_card(group_label, [])


CSS_BLOCK = """
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --text-color: #f8fafc;
            --border-color: #334155;
            --primary-color: #60a5fa;
        }

        html, body {
            margin: 0;
            padding: 0;
            background-color: var(--bg-color);
            color: var(--text-color);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            box-sizing: border-box;
        }

        .dashboard-container {
            display: flex;
            flex-direction: column;
            padding: 16px 20px 40px;
            box-sizing: border-box;
        }

        header {
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 10px;
            margin-bottom: 16px;
            gap: 16px;
            flex-wrap: wrap;
        }

        .header-text h1 { margin: 0; font-size: 1.4rem; color: #f1f5f9; }
        .header-text p.subtitle { margin: 4px 0 0; color: #94a3b8; font-size: 0.85rem; }

        .header-actions {
            display: flex;
            gap: 8px;
            flex-shrink: 0;
        }

        .download-btn {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: var(--card-bg);
            color: var(--primary-color);
            border: 1px solid var(--primary-color);
            padding: 7px 14px;
            border-radius: 6px;
            text-decoration: none;
            font-size: 0.8rem;
            font-weight: 600;
            font-family: inherit;
            cursor: pointer;
            white-space: nowrap;
        }
        .download-btn:hover { background: var(--primary-color); color: #0f172a; }

        .sections-stack {
            display: flex;
            flex-direction: column;
            gap: 22px;
        }

        .section {
            background: var(--card-bg);
            border-radius: 8px;
            padding: 14px 16px;
            border: 1px solid var(--border-color);
        }

        .section-title {
            font-size: 1.1rem;
            margin: 0 0 12px 0;
            color: var(--primary-color);
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 6px;
        }

        .panel { margin-bottom: 18px; }
        .panel:last-child { margin-bottom: 0; }

        .genus-subgroup { margin-bottom: 14px; }
        .genus-subgroup:last-child { margin-bottom: 0; }

        .genus-subheading {
            font-size: 0.78rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: #94a3b8;
            margin-bottom: 6px;
            padding-bottom: 4px;
            border-bottom: 1px dashed var(--border-color);
        }

        .panel-label {
            font-size: 0.95rem;
            font-weight: 600;
            color: #cbd5e1;
            margin-bottom: 8px;
        }

        .group-row {
            display: flex;
            flex-wrap: wrap;
            justify-content: center;
            gap: 10px;
        }

        /* A row with too few cards to naturally fill the available width
           would otherwise sit at the fixed 260-420px base card size,
           leaving a large empty gap and making the dashboard look sparse
           for small comparisons (a genus with only 2-3 species). For rows
           of 2..MAX_FILL_COLUMNS cards (MAX_FILL_COLUMNS is defined in
           Python; --fill-n is set per-row from the actual card count),
           each card's width is computed with calc() so that exactly n of
           them tile the row edge to edge, accounting for the 10px gaps --
           a 2-card row splits the row evenly in half, a 3-card row into
           thirds, and so on, with zero leftover space regardless of the
           row's actual pixel width. A single-card row (a single genus with
           only one match) is deliberately NOT stretched edge to edge --
           that looked odd (one image blown up to the full row width) -- it
           instead sits at the normal 260-420px card size, centered by the
           .group-row justify-content above. Rows above MAX_FILL_COLUMNS
           cards also keep the base sizing and simply wrap onto more lines,
           since stretching e.g. 8 cards edge-to-edge each would make them
           too large for comfortable side-by-side comparison. */
        .group-row.fill-row > .group-card {
            flex: 1 1 calc((100% - (var(--fill-n) - 1) * 10px) / var(--fill-n));
            max-width: calc((100% - (var(--fill-n) - 1) * 10px) / var(--fill-n));
        }
        .group-row.fill-row > .group-card img {
            max-height: 650px;
        }

        .group-card {
            box-sizing: border-box;
            flex: 1 1 260px;
            max-width: 420px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            overflow: hidden;
            background: #020617;
            display: flex;
            flex-direction: column;
        }

        .group-card-label {
            font-size: 0.8rem;
            font-weight: 700;
            color: var(--primary-color);
            background: #0b1220;
            padding: 4px 8px;
            border-bottom: 1px solid var(--border-color);
        }

        .group-card a { display: flex; justify-content: center; align-items: center; padding: 4px; }

        .group-card img {
            max-width: 100%;
            max-height: 320px;
            width: auto;
            height: auto;
            object-fit: contain;
            cursor: zoom-in;
        }

        .card-caption {
            padding: 2px 8px;
            font-size: 0.68rem;
            color: #94a3b8;
            background: #0f172a;
            border-top: 1px solid var(--border-color);
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }

        .missing-card {
            justify-content: center;
            align-items: center;
            padding: 14px 8px;
            border: 1px dashed #f87171;
            background: #450a0a;
        }
        .missing-msg { color: #fca5a5; font-size: 0.8rem; padding-top: 6px; }

        .reference-tree-section {
            border: 1px solid var(--primary-color);
        }
        .reference-tree-card {
            display: flex;
            justify-content: center;
            background: #020617;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 10px;
        }
        .reference-tree-card img {
            max-width: 100%;
            max-height: 480px;
            width: auto;
            height: auto;
            object-fit: contain;
            cursor: zoom-in;
        }

        .kpi-row {
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            margin-bottom: 18px;
        }
        .kpi-card {
            flex: 1 1 130px;
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 12px 10px;
            text-align: center;
        }
        .kpi-value {
            font-size: 1.6rem;
            font-weight: 800;
            color: var(--primary-color);
            line-height: 1.1;
        }
        .kpi-label {
            font-size: 0.72rem;
            color: #94a3b8;
            margin-top: 4px;
        }

        .arch-map-section { border: 1px solid var(--border-color); }
        .arch-map-svg {
            width: 100%;
            height: auto;
            background: #ffffff;
            border: 1px solid var(--border-color);
            border-radius: 6px;
        }

        .lightbox-overlay {
            display: none;
            position: fixed;
            top: 0; left: 0;
            width: 100%; height: 100%;
            background: rgba(2, 6, 23, 0.92);
            z-index: 1000;
            justify-content: center;
            align-items: center;
            flex-direction: column;
            padding: 24px;
            box-sizing: border-box;
        }
        .lightbox-overlay.active { display: flex; }

        /* Scrollable viewport: the image is NO LONGER permanently shrunk to
           fit. At "Fit" zoom it scales down as before, but at any zoom above
           that the viewport scrolls in BOTH directions (overflow: auto), so
           every part of a large figure is reachable at full detail instead
           of being squashed into 78vh. */
        .lightbox-viewport {
            width: 92vw;
            height: 76vh;
            overflow: auto;
            background: #fff;
            border: 2px solid var(--border-color);
            border-radius: 6px;
            display: flex;
            align-items: center;
            justify-content: center;
        }
        /* Once zoomed past fit, anchor top-left so scrolling starts at the
           image's origin rather than centering it and clipping both edges. */
        .lightbox-viewport.zoomed {
            align-items: flex-start;
            justify-content: flex-start;
        }

        .lightbox-content {
            display: block;
            max-width: 100%;
            max-height: 100%;
            object-fit: contain;
        }
        .lightbox-viewport.zoomed .lightbox-content {
            max-width: none;
            max-height: none;
            height: auto;
        }

        .lightbox-caption { margin-top: 10px; color: #94a3b8; font-size: 0.85rem; }

        .lightbox-toolbar {
            margin-top: 14px;
            display: flex;
            gap: 10px;
            align-items: center;
            flex-wrap: wrap;
            justify-content: center;
        }

        .lightbox-zoom-level {
            color: #cbd5e1;
            font-size: 0.85rem;
            min-width: 72px;
            text-align: center;
            font-variant-numeric: tabular-nums;
        }

        /* Inline SVG (Promoter Architecture Map) is generated without a
           src, unlike every other panel's base64 <img> -- clicking
           anywhere on the figure rasterizes it on the fly and opens the
           same lightbox everything else uses. See wrap_arch_map_figure. */
        .arch-map-figure { position: relative; cursor: zoom-in; }

        /* The map always scales to fit its container width -- no
           horizontal scrolling, so the whole map is visible at once. */
        .arch-map-scroll { overflow: visible; }

        .lightbox-btn {
            background: var(--primary-color);
            color: #0f172a;
            padding: 9px 20px;
            border-radius: 6px;
            text-decoration: none;
            font-weight: 700;
            font-size: 0.9rem;
            cursor: pointer;
            border: none;
            font-family: inherit;
        }

        .lightbox-close-btn {
            background: transparent;
            color: #f8fafc;
            border: 1px solid var(--border-color);
        }
"""

LIGHTBOX_HTML_AND_JS = """
    <div class="lightbox-overlay" id="lightbox-overlay" onclick="if (event.target === this) closeLightbox();">
        <div class="lightbox-viewport" id="lightbox-viewport">
            <img class="lightbox-content" id="lightbox-img" src="" alt="Full size view"/>
        </div>
        <div class="lightbox-caption" id="lightbox-caption"></div>
        <div class="lightbox-toolbar">
            <button class="lightbox-btn lightbox-close-btn" onclick="zoomLightbox(-1)">&minus;</button>
            <span class="lightbox-zoom-level" id="lightbox-zoom-level">Fit</span>
            <button class="lightbox-btn lightbox-close-btn" onclick="zoomLightbox(1)">+</button>
            <button class="lightbox-btn lightbox-close-btn" onclick="setLightboxZoom(null)">Fit</button>
            <button class="lightbox-btn lightbox-close-btn" onclick="setLightboxZoom(1)">100%</button>
            <a class="lightbox-btn" id="lightbox-download" href="" download="image.png">&#8681; Download image</a>
            <button class="lightbox-btn lightbox-close-btn" onclick="closeLightbox()">&times; Close</button>
        </div>
    </div>

    <script>
        // null = "fit to viewport" (default). A number = explicit scale
        // factor applied to the image's NATURAL pixel size, at which point
        // the viewport scrolls in both directions to reach any part of it.
        var _lbZoom = null;
        var _LB_STEPS = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8];

        function _applyLightboxZoom() {
            var img = document.getElementById('lightbox-img');
            var viewport = document.getElementById('lightbox-viewport');
            var label = document.getElementById('lightbox-zoom-level');
            if (_lbZoom === null) {
                viewport.classList.remove('zoomed');
                img.style.width = '';
                label.textContent = 'Fit';
            } else {
                viewport.classList.add('zoomed');
                // naturalWidth is 0 until the image decodes; _applyLightboxZoom
                // is re-called from the onload handler below for that case.
                if (img.naturalWidth) {
                    img.style.width = (img.naturalWidth * _lbZoom) + 'px';
                }
                label.textContent = Math.round(_lbZoom * 100) + '%';
            }
        }

        function setLightboxZoom(z) {
            _lbZoom = z;
            _applyLightboxZoom();
        }

        function zoomLightbox(direction) {
            var current = _lbZoom;
            if (current === null) {
                // Stepping away from "Fit" starts from the actual on-screen
                // scale, so "+" visibly enlarges rather than jumping to a
                // smaller fixed step than the current fitted size.
                var img = document.getElementById('lightbox-img');
                current = (img.naturalWidth && img.clientWidth)
                    ? (img.clientWidth / img.naturalWidth) : 1;
            }
            var i = 0;
            while (i < _LB_STEPS.length && _LB_STEPS[i] <= current + 0.001) { i++; }
            var next = direction > 0 ? _LB_STEPS[Math.min(i, _LB_STEPS.length - 1)]
                                     : _LB_STEPS[Math.max(0, i - 2)];
            setLightboxZoom(next);
        }

        function openLightbox(src, filename) {
            var img = document.getElementById('lightbox-img');
            img.onload = _applyLightboxZoom;   // natural size known only after decode
            img.src = src;
            document.getElementById('lightbox-caption').textContent = filename;
            var dl = document.getElementById('lightbox-download');
            dl.href = src;
            dl.download = filename;
            setLightboxZoom(null);             // always open at Fit
            document.getElementById('lightbox-overlay').classList.add('active');
        }

        function closeLightbox() {
            document.getElementById('lightbox-overlay').classList.remove('active');
        }

        document.addEventListener('keydown', function (e) {
            if (!document.getElementById('lightbox-overlay').classList.contains('active')) return;
            if (e.key === 'Escape') closeLightbox();
            if (e.key === '+' || e.key === '=') zoomLightbox(1);
            if (e.key === '-' || e.key === '_') zoomLightbox(-1);
            if (e.key === '0') setLightboxZoom(null);
        });

        // Ctrl/Cmd + wheel zooms; a plain wheel still scrolls the viewport
        // normally, which is what you want for panning a zoomed figure.
        document.getElementById('lightbox-viewport').addEventListener('wheel', function (e) {
            if (!(e.ctrlKey || e.metaKey)) return;
            e.preventDefault();
            zoomLightbox(e.deltaY < 0 ? 1 : -1);
        }, {passive: false});

        // ---- Inline-SVG (Promoter Architecture Map) zoom helper ----
        // Every other panel is a base64 <img>, so the lightbox's own
        // download link covers it for free. The architecture map is
        // generated as INLINE SVG markup instead, which has no src to hand
        // to that link -- so clicking it rasterizes the live SVG on the
        // fly and opens it in the same zoomable lightbox as everything else.

        function _svgElementFrom(btn) {
            return btn.closest('.arch-map-figure').querySelector('svg');
        }

        // A clone with explicit width/height from its viewBox. The SVG is
        // authored with only a viewBox (so it scales responsively in-page),
        // but an Image() used for canvas drawing needs real intrinsic
        // dimensions or it renders at zero size.
        function _svgCloneWithSize(svgEl) {
            var vb = (svgEl.getAttribute('viewBox') || '').split(/[\\s,]+/).map(Number);
            var w = (vb.length === 4 && vb[2]) ? vb[2] : (svgEl.clientWidth || 780);
            var h = (vb.length === 4 && vb[3]) ? vb[3] : (svgEl.clientHeight || 300);
            var clone = svgEl.cloneNode(true);
            clone.setAttribute('width', w);
            clone.setAttribute('height', h);
            clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
            return {clone: clone, width: w, height: h};
        }

        // PNG via canvas. Safe here because this SVG is fully self-contained
        // (literal hex fills, no CSS variables, no external refs, no embedded
        // raster) -- so the canvas never gets tainted and nothing renders
        // as blank. Rendered at 3x so the lightbox's zoom has real detail.
        function _archMapPngDataUrl(btn, scale, callback) {
            var sized = _svgCloneWithSize(_svgElementFrom(btn));
            var markup = new XMLSerializer().serializeToString(sized.clone);
            var svgUrl = 'data:image/svg+xml;base64,'
                       + btoa(unescape(encodeURIComponent(markup)));
            var img = new Image();
            img.onload = function () {
                var canvas = document.createElement('canvas');
                canvas.width = sized.width * scale;
                canvas.height = sized.height * scale;
                var ctx = canvas.getContext('2d');
                ctx.fillStyle = '#ffffff';   // SVG has no background of its own
                ctx.fillRect(0, 0, canvas.width, canvas.height);
                ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
                callback(canvas.toDataURL('image/png'));
            };
            img.onerror = function () {
                alert('Could not open the architecture map at full size.');
            };
            img.src = svgUrl;
        }

        // Opens the architecture map in the same zoomable/scrollable
        // lightbox as every other panel, so it gets identical zoom + the
        // lightbox's own download link.
        function openArchMapLightbox(btn, filename) {
            _archMapPngDataUrl(btn, 3, function (dataUrl) {
                openLightbox(dataUrl, filename);
            });
        }
    </script>
"""


DOWNLOAD_BUTTON_HTML = (
    '<div class="header-actions">'
    '<button class="download-btn" onclick="downloadSpeciesDashboard()">'
    '&#8681; Download species-wise dashboard</button>'
    '</div>'
)

DOWNLOAD_SCRIPT_TEMPLATE = """
    <script>
        function downloadSpeciesDashboard() {{
            var blob = new Blob([SPECIES_DASHBOARD_HTML], {{type: 'text/html'}});
            var url = URL.createObjectURL(blob);
            var a = document.createElement('a');
            a.href = url;
            a.download = {species_filename_json};
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            setTimeout(function() {{ URL.revokeObjectURL(url); }}, 1000);
        }}
        const SPECIES_DASHBOARD_HTML = {species_html_json};
    </script>
"""


def render_sections_html(sections, groups, group_by_genus=False):
    """Renders just the section/panel/group-card markup for one dashboard
    variant (genus or species) -- shared by both the visible page and the
    embedded standalone document, so there's exactly one place that turns
    `sections` + `groups` into HTML. Returns (html_fragment, n_found, n_missing).

    group_by_genus=True (used only for the species-wise dashboard) splits
    each panel's row of species cards into one sub-row PER GENUS, with a
    small genus sub-heading above each -- so e.g. every Plasmodium species
    is visually grouped together and separated from every Toxoplasma
    species, instead of all species from every genus being mixed into one
    flat row. The genus for each species is derived the same way
    find_images_for_group does (first underscore-separated token of the
    sanitized species id), so it lines up with how patterns are matched."""
    # Rows of 2..this many cards get the "fill-row" treatment (see the
    # .group-row.fill-row CSS rule) so they tile the full row width with no
    # empty gap, instead of sitting at the fixed 260-420px base card size.
    # A single-card row is deliberately excluded -- it centers at the base
    # card size instead (see the CSS comment) rather than being stretched
    # edge to edge. This applies uniformly to genus rows AND to each
    # genus's species sub-row in the species-wise dashboard -- both go
    # through render_group_row below.
    MAX_FILL_COLUMNS = 4

    normalized_groups = [g if isinstance(g, tuple) else (g, g) for g in groups]
    n_found, n_missing = 0, 0
    parts = []

    genus_buckets = None
    if group_by_genus:
        genus_buckets = {}
        for sid, label in normalized_groups:
            genus = str(sid).split("_")[0]
            genus_buckets.setdefault(genus, []).append((sid, label))

    def render_group_row(group_list):
        nonlocal n_found, n_missing
        n = len(group_list)
        if 2 <= n <= MAX_FILL_COLUMNS:
            row_class = "group-row fill-row"
            row_style = f' style="--fill-n:{n};"'
        else:
            row_class = "group-row"
            row_style = ""
        row_parts = [f"""
                <div class="{row_class}"{row_style}>"""]
        for group_id, group_label in group_list:
            matches = find_images_for_group(
                group_id, panel["patterns"], OUTPUT_DIR, SHOW_ALL_SPECIES_PER_GENUS
            )
            if matches:
                n_found += 1
            else:
                n_missing += 1
            row_parts.append(render_group_card(group_label, matches))
        row_parts.append("""
                </div>""")
        return "".join(row_parts)

    for section in sections:
        parts.append(f"""
        <div class="section">
            <h2 class="section-title">{section["title"]}</h2>""")
        for panel in section["panels"]:
            parts.append(f"""
            <div class="panel">
                <div class="panel-label">{panel["label"]}</div>""")

            if group_by_genus:
                for genus, genus_species in genus_buckets.items():
                    parts.append(f"""
                <div class="genus-subgroup">
                    <div class="genus-subheading">{genus}</div>{render_group_row(genus_species)}
                </div>""")
            else:
                parts.append(render_group_row(normalized_groups))

            parts.append("""
            </div>""")
        parts.append("""
        </div>""")

    return "".join(parts), n_found, n_missing


REFERENCE_SPECIES_TREE_PATTERNS = [
    "reference_species_tree/orthofinder_species_tree.png",
]

# The 03c tree-building script now writes ONE tree per genus (when it
# succeeds for that genus) under this layout, in addition to attempting the
# single "all genera combined" tree above. These per-genus patterns are
# tried per detected genus so a genus like Toxoplasma -- which built fine --
# still shows its tree even when the combined, all-genera tree failed (e.g.
# too few orthogroups, as happened for Babesia in this run).
REFERENCE_SPECIES_TREE_PATTERNS_BY_GENUS = [
    "reference_species_tree/by_genus/{genus}/{genus}_orthofinder_species_tree.png",
]


# ============================================================================
# NEW: Summary KPI cards + Promoter Architecture Map
#
# Both are built from REAL parsed data (MEME's own meme.xml, TOMTOM's own
# tomtom.tsv, the actual promoter FASTA, and the CpG detailed CSVs) -- not
# guessed from filenames the way most panels above are, since these need
# actual VALUES out of the files, not just "does this file exist".
#
# WHY MOTIF <-> TOMTOM MATCHING IS DONE BY ORDINAL POSITION, NOT BY STRING
# ----------------------------------------------------------------------------
# meme.xml's <motif id="motif_1" ...> uses an internal ID that does NOT
# necessarily match tomtom.tsv's Query_ID column verbatim (MEME versions
# differ on whether that's "1", "motif_1", "MEME-1", etc.), and I don't have
# a real tomtom.tsv sample from this pipeline to confirm which. But TOMTOM
# is always run against the SAME meme.xml, motif-by-motif, in the SAME
# order MEME assigned them -- so the i-th <motif> element in meme.xml
# corresponds to the i-th DISTINCT Query_ID that first appears in
# tomtom.tsv, regardless of the exact string format either one uses. This
# was verified against a schema-accurate fixture before being wired in
# here (see the conversation this was built in) -- if your real tomtom.tsv
# ever has motifs in a different order than meme.xml (shouldn't happen in
# normal MEME Suite usage, but flagging it), this pairing would be wrong.
# ============================================================================

def parse_meme_xml(meme_xml_path):
    """
    Returns (seq_info, motifs):
      seq_info: {internal_seq_id: {"name": real_fasta_header, "length": int}}
      motifs:   {internal_motif_id: {"name": consensus, "alt": alt_name,
                                       "width": int, "sites": [site, ...]}}
      each site: {"seq_id": internal_seq_id, "position": int (1-based,
                  per the MEME XML DTD), "strand": "plus"/"minus",
                  "pvalue": float}
    Returns ({}, {}) if the file doesn't exist or fails to parse -- never
    raises, since a missing/malformed meme.xml for one organism shouldn't
    break the whole dashboard.
    """
    if not os.path.exists(meme_xml_path):
        return {}, {}
    try:
        root = ET.parse(meme_xml_path).getroot()
    except ET.ParseError:
        return {}, {}

    seq_info = {}
    training_set = root.find("training_set")
    if training_set is not None:
        for seq in training_set.findall("sequence"):
            seq_info[seq.get("id")] = {
                "name": seq.get("name"),
                "length": int(seq.get("length", 0)),
            }

    motifs = {}
    motifs_el = root.find("motifs")
    if motifs_el is not None:
        for motif in motifs_el.findall("motif"):
            sites = []
            contributing = motif.find("contributing_sites")
            if contributing is not None:
                for site in contributing.findall("contributing_site"):
                    try:
                        sites.append({
                            "seq_id": site.get("sequence_id"),
                            "position": int(site.get("position", 0)),
                            "strand": site.get("strand"),
                            "pvalue": float(site.get("pvalue", 1.0)),
                        })
                    except (TypeError, ValueError):
                        continue
            motifs[motif.get("id")] = {
                "name": motif.get("name"),
                "alt": motif.get("alt"),
                "width": int(motif.get("width", 0) or 0),
                "sites": sites,
            }
    return seq_info, motifs


def strip_motif_file_extension(name):
    """Drop a trailing '.bml' (case-insensitive) from a motif/target name,
    e.g. 'PF11_0442_RC.bml' -> 'PF11_0442_RC'. Other names are unchanged."""
    return re.sub(r"\.bml$", "", str(name), flags=re.IGNORECASE)


def parse_tomtom_tsv(tomtom_tsv_path):
    """
    Returns an ORDERED list of (query_id_raw, best_target_id, best_qvalue)
    -- one entry per distinct Query_ID, in the order each first appears,
    keeping only the best (lowest q-value) match per query when a motif
    hits multiple targets. Skips MEME's trailing "#"-prefixed comment
    lines and blank lines, both of which real tomtom.tsv output has at
    the end of the file.
    """
    if not os.path.exists(tomtom_tsv_path):
        return []
    with open(tomtom_tsv_path, newline="") as f:
        lines = [ln for ln in f if ln.strip() and not ln.startswith("#")]
    if not lines:
        return []

    best_by_query = {}
    order = []
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        qid = row.get("Query_ID")
        if qid is None:
            continue
        try:
            qval = float(row.get("q-value", "1"))
        except (TypeError, ValueError):
            continue
        if qid not in best_by_query or qval < best_by_query[qid][1]:
            best_by_query[qid] = (strip_motif_file_extension(row.get("Target_ID", "?")), qval)
        if qid not in order:
            order.append(qid)

    return [(qid, best_by_query[qid][0], best_by_query[qid][1]) for qid in order]


def summarize_group_motifs(meme_dir, tomtom_dir, group_label):
    """
    For one organism/genus: parses its meme.xml, matches each motif to its
    best TOMTOM hit (across every reference database subfolder present --
    e.g. JASPAR2026_CORE and UniProbe -- keeping the single best q-value
    seen for that motif across all of them), and collapses each motif's
    scattered occurrences into ONE representative bar: its MEDIAN
    TSS-relative start position (median, not mean, so one extreme outlier
    occurrence doesn't drag the bar to a misleading position) and its
    width (fixed per motif, from meme.xml).

    Returns a list of dicts: {"label", "start_tss", "width", "n_sites",
    "matched"} -- one per motif that had at least one valid occurrence.
    Empty list if meme.xml is missing/empty for this group.
    """
    meme_xml_path = os.path.join(meme_dir, group_label, "meme.xml")
    seq_info, motifs = parse_meme_xml(meme_xml_path)
    if not motifs:
        return []

    motif_ids_in_order = list(motifs.keys())

    # ---- resolve each TOMTOM Query_ID to the meme.xml motif it belongs to ----
    # tomtom.tsv only contains motifs that had >= 1 hit, so counting rows
    # ("i-th query == i-th motif") shifts every name onto the wrong motif as
    # soon as one motif has no hit. Match on the Query_ID text instead
    # (motif id / name / alt / "MEME-3" / "motif_3" / "3-CONSENSUS" ...).
    def _norm(x):
        return re.sub(r"[^a-z0-9]", "", str(x or "").lower())

    key_to_ordinal = {}
    for ordinal, motif_id in enumerate(motif_ids_in_order):
        mm = motifs[motif_id]
        for k in (motif_id, mm.get("name"), mm.get("alt")):
            nk = _norm(k)
            if nk:
                key_to_ordinal.setdefault(nk, ordinal)

    def resolve_query(qid, n_queries):
        nk = _norm(qid)
        if nk in key_to_ordinal:
            return key_to_ordinal[nk]
        for pat in (r"^(?:meme|motif)?[-_ ]?(\d+)$", r"^(\d+)[-_ ]", r"(\d+)$"):
            mt = re.search(pat, str(qid).strip(), flags=re.I)
            if mt and 1 <= int(mt.group(1)) <= len(motif_ids_in_order):
                return int(mt.group(1)) - 1
        return None

    best_match_by_ordinal = {}
    group_tomtom_dir = os.path.join(tomtom_dir, group_label)
    pairing_modes = []
    if os.path.isdir(group_tomtom_dir):
        for db_name in sorted(os.listdir(group_tomtom_dir)):
            tsv_path = os.path.join(group_tomtom_dir, db_name, "tomtom.tsv")
            queries = parse_tomtom_tsv(tsv_path)
            if not queries:
                continue
            # Preferred: pair by Query_ID text -- only used when EVERY query in
            # this file resolves to a distinct motif, so it can never be
            # partly right and partly wrong.
            resolved = [resolve_query(qid, len(queries)) for qid, _t, _q in queries]
            by_id = (all(r is not None for r in resolved)
                     and len(set(resolved)) == len(resolved))
            # Otherwise: EXACTLY the original behaviour (i-th distinct
            # Query_ID in the file == i-th motif in meme.xml), so this never
            # names fewer motifs than the previous version of the script did.
            pairing_modes.append(f"{db_name}:{'query-id' if by_id else 'positional'}")
            for pos, (qid, target, qval) in enumerate(queries):
                ordinal = resolved[pos] if by_id else pos
                if ordinal is None or ordinal >= len(motif_ids_in_order):
                    continue
                if ordinal not in best_match_by_ordinal or qval < best_match_by_ordinal[ordinal][1]:
                    best_match_by_ordinal[ordinal] = (target, qval)

    results = []
    for ordinal, motif_id in enumerate(motif_ids_in_order):
        m = motifs[motif_id]
        starts_tss = []
        for site in m["sites"]:
            seq = seq_info.get(site["seq_id"])
            if not seq or not seq["length"]:
                continue
            # MEME positions are 1-BASED (first base = 1), so the last base of
            # the promoter (position == length) must map to -1, not 0:
            # TSS-relative = position - length - 1. (Same convention as 07f/07i.)
            starts_tss.append(site["position"] - seq["length"] - 1)
        if not starts_tss:
            continue
        starts_tss.sort()
        median_start = starts_tss[len(starts_tss) // 2]

        match = best_match_by_ordinal.get(ordinal)
        if match and match[1] <= TOMTOM_QVALUE_CUTOFF:
            label, matched = match[0], True
        else:
            label, matched = f"Motif {ordinal + 1}", False

        results.append({
            "label": label,
            "start_tss": median_start,
            "width": m["width"],
            "n_sites": len(starts_tss),
            "matched": matched,
            "meme_motif": f"Motif {ordinal + 1}",
            "qvalue": match[1] if match else None,
        })
    print(f"[architecture map] {group_label}: {len(results)} motif(s), "
          f"{sum(1 for r in results if r['matched'])} matched to TOMTOM (q <= {TOMTOM_QVALUE_CUTOFF}); "
          f"pairing [{', '.join(pairing_modes) or 'no tomtom.tsv found'}]; best q per motif: "
          + ", ".join(f"M{o + 1}={best_match_by_ordinal[o][1]:.3g}" for o in sorted(best_match_by_ordinal)))
    return results


def count_fasta_records(fasta_path):
    if not os.path.exists(fasta_path):
        return 0
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith(">"))


def count_csv_data_rows(csv_path):
    if not os.path.exists(csv_path):
        return 0
    with open(csv_path, newline="") as f:
        n = sum(1 for _ in csv.reader(f))
    return max(0, n - 1)  # minus header row


def count_total_meme_motifs(meme_output_dir):
    """Sums <motif> counts across every organism's meme.xml under
    meme_output_dir -- iterates actual folders present rather than a
    precomputed species list, so it's correct even if a species is
    missing its own MEME run."""
    if not os.path.isdir(meme_output_dir):
        return 0
    total = 0
    for entry in sorted(os.listdir(meme_output_dir)):
        _, motifs = parse_meme_xml(os.path.join(meme_output_dir, entry, "meme.xml"))
        total += len(motifs)
    return total


def count_total_cpg_islands(cpg_by_organism_dir):
    """Sums CpG island counts from each organism's *_cpg_islands_detailed.csv
    (one row per island) -- NOT also summing the _By_Genus versions, since
    those are the same islands re-pooled by genus and would double-count."""
    if not os.path.isdir(cpg_by_organism_dir):
        return 0
    total = 0
    for entry in sorted(os.listdir(cpg_by_organism_dir)):
        csv_path = os.path.join(cpg_by_organism_dir, entry, f"{entry}_cpg_islands_detailed.csv")
        total += count_csv_data_rows(csv_path)
    return total


def compute_kpis():
    """Returns an ordered list of (value, label) for the KPI card row.
    Every value comes from actually parsing/counting real output files --
    see each counter function's docstring for exactly which file(s)."""
    n_species = len(detect_species(LABELS_PATH))
    n_promoters = count_fasta_records(PROMOTER_FASTA_PATH)
    n_tf_motifs = count_total_meme_motifs(MEME_OUTPUT_DIR_ORGANISM)
    n_cpg_islands = count_total_cpg_islands(CPG_DETAILED_DIR_ORGANISM)
    return [
        (n_species, "Species"),
        (n_promoters, "Promoters"),
        (n_tf_motifs, "TF Motifs (MEME-discovered)"),
        (n_cpg_islands, "CpG Islands"),
    ]


def render_kpi_cards_html():
    kpis = compute_kpis()
    cards = "".join(
        f'<div class="kpi-card"><div class="kpi-value">{value}</div>'
        f'<div class="kpi-label">{label}</div></div>'
        for value, label in kpis
    )
    return f'<div class="kpi-row">{cards}</div>'


CORE_PROMOTER_RESULTS_DIR_ORGANISM = os.path.join(OUTPUT_DIR, "Core_Promoter_Results_By_Organism")
# 06b actually writes GENUS-pooled results (subfolders "Plasmodium", "Toxoplasma", ...)
# into the "_By_Organism" folder above. If the stage is ever fixed/renamed to emit
# Core_Promoter_Results_By_Genus/, that folder is preferred for every genus-level
# lookup; the legacy folder is used only when it doesn't exist. The legacy folder is
# still tried FIRST for exact per-organism subfolders, and anything resolved from a
# genus folder is flagged genus_pooled so the caption says so.
CORE_PROMOTER_RESULTS_DIR_GENUS = os.path.join(OUTPUT_DIR, "Core_Promoter_Results_By_Genus")


def _core_promoter_genus_dir():
    if os.path.isdir(CORE_PROMOTER_RESULTS_DIR_GENUS):
        return CORE_PROMOTER_RESULTS_DIR_GENUS
    return CORE_PROMOTER_RESULTS_DIR_ORGANISM

# Fixed, recognizable colors for named core promoter elements -- these are
# a small, well-known fixed set (unlike arbitrary MEME-discovered motifs,
# which get auto-cycled colors instead), matched by substring so "TATA
# box" and "Initiator (Inr)" both hit their entry regardless of the exact
# label text a given run uses.
CORE_ELEMENT_COLORS = {
    "tata": "#f97316",   # orange
    "inr":  "#8b5cf6",   # purple
    "bre":  "#14b8a6",   # teal
    "dpe":  "#ec4899",   # pink
    "mte":  "#84cc16",   # lime
    "caat": "#eab308",   # yellow
}
CORE_ELEMENT_FALLBACK_COLOR = "#94a3b8"

# Short descriptions for the same fixed set of named core elements, shown
# under each element's box in the architecture map -- matches the textbook
# core-promoter diagram convention (e.g. BRE -> "TFIIB recognition
# element"). Matched by the same substring rule as CORE_ELEMENT_COLORS.
# MEME-discovered motifs have no entry here and simply get no description
# line under their box, since there's no fixed descriptive text for an
# arbitrary discovered motif the way there is for a named core element.
CORE_ELEMENT_DESCRIPTIONS = {
    "tata": "TATA box",
    "inr":  "Initiator",
    "bre":  "TFIIB recognition element",
    "dpe":  "Downstream promoter element",
    "mte":  "Motif ten element",
    "caat": "CCAAT box",
}


def color_for_core_element(label):
    key = str(label).lower()
    for token, color in CORE_ELEMENT_COLORS.items():
        if token in key:
            return color
    return CORE_ELEMENT_FALLBACK_COLOR


def description_for_core_element(label):
    key = str(label).lower()
    for token, desc in CORE_ELEMENT_DESCRIPTIONS.items():
        if token in key:
            return desc
    return ""


def _contrast_text_color(hex_color):
    """Picks black or white text so a label stays readable against any of
    this map's fill colors (fixed core-element colors, auto-cycled motif
    colors, and the grey 'unmatched' color) without hand-tuning per color."""
    hex_color = str(hex_color).lstrip("#")
    try:
        r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    except (ValueError, IndexError):
        return "#000000"
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "#000000" if luminance > 150 else "#ffffff"


def load_core_promoter_hits(csv_path):
    """Reads a *_core_promoter_individual_hits.csv (confirmed real columns:
    seq_id, motif, strand, start_pos_0based, rel_tss_pos_bp,
    matched_sequence). Returns {motif_name: [(rel_tss_pos_bp, width), ...]}.
    Missing file -> {} (not every organism necessarily has this stage's
    output)."""
    if not os.path.exists(csv_path):
        return {}
    hits = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            motif = row.get("motif")
            if not motif:
                continue
            try:
                pos = int(float(row.get("rel_tss_pos_bp")))
            except (TypeError, ValueError):
                continue
            width = len(row.get("matched_sequence") or "") or 6
            hits.setdefault(motif, []).append((pos, width))
    return hits


def _collapse_core_hits(hits):
    """Shared collapse logic: median position + median width per element,
    same "median not mean" reasoning as summarize_group_motifs()."""
    results = []
    for motif, entries in hits.items():
        positions = sorted(p for p, _ in entries)
        widths = sorted(w for _, w in entries)
        results.append({
            "label": motif,
            "start_tss": positions[len(positions) // 2],
            "width": widths[len(widths) // 2],
            "n_sites": len(entries),
            "matched": True,          # named core elements are never "unmatched"
            "is_core_element": True,
        })
    return results


def summarize_core_promoter_elements_for_group(group_label, is_genus):
    """
    Returns collapsed core-element bars (BRE/TATA/Inr/DPE/MTE/CAAT) for one
    group, read from Core_Promoter_Results_By_Organism's
    *_core_promoter_individual_hits.csv.

    IMPORTANT NAMING GOTCHA this function has to work around: despite its
    "_By_Organism" name, Core_Promoter_Results_By_Organism/ is keyed by
    GENUS (its subfolders are "Plasmodium", "Toxoplasma", ... not
    "Plasmodium_falciparum_3D7"), because 06b's organism-wise stage groups
    by genus. That's confirmed by the actual outputs/ listing, and it's
    the same situation the image panels in SPECIES_SECTIONS already handle
    with their {genus} fallback pattern.

    is_genus=True: reads the genus folder directly (the startswith() scan
    also catches any truly per-organism folders if that stage is ever
    fixed to emit them).
    is_genus=False (species-wise dashboard): tries that species' own folder
    FIRST -- so this automatically starts showing real per-species
    positions if 06b is ever made truly per-organism -- then falls back to
    its parent genus's folder. Without that fallback this returned {} for
    every species, which is why no BRE/TATA/Inr bars appeared on the
    species-wise architecture map at all.

    NOTE ON WHAT THE FALLBACK ACTUALLY MEANS: when the genus fallback is
    used, every species in a genus shows the SAME genus-pooled element
    positions -- it is not independent per-species evidence. The rendered
    caption says so explicitly rather than letting identical bars across
    congeneric species read as independent confirmation.
    """
    genus_dir = _core_promoter_genus_dir()
    if not (os.path.isdir(CORE_PROMOTER_RESULTS_DIR_ORGANISM) or os.path.isdir(genus_dir)):
        return []

    def read_folder(entry, base_dir=None):
        base_dir = base_dir or CORE_PROMOTER_RESULTS_DIR_ORGANISM
        csv_path = os.path.join(base_dir, entry,
                                 f"{entry}_core_promoter_individual_hits.csv")
        return load_core_promoter_hits(csv_path)

    if not is_genus:
        own = read_folder(group_label)
        if own:
            return _collapse_core_hits(own)
        genus = str(group_label).split("_")[0]
        pooled = read_folder(genus, genus_dir)
        if pooled:
            results = _collapse_core_hits(pooled)
            for r in results:
                r["genus_pooled"] = True
            return results
        return []

    combined = {}
    for entry in sorted(os.listdir(genus_dir)):
        if not entry.startswith(group_label):
            continue
        for motif, entries in read_folder(entry, genus_dir).items():
            combined.setdefault(motif, []).extend(entries)
    return _collapse_core_hits(combined)


# Palette for TOMTOM-matched motifs -- deliberately different from the
# core-element colours (teal/yellow/pink/purple/lime/orange) so a matched
# TF never looks like BREu, CAAT, DPE, Inr, MTE or TATA.
ARCH_MAP_COLORS = ["#2563eb", "#dc2626", "#0891b2", "#78350f", "#16a34a",
                   "#1e3a8a", "#9f1239", "#4d7c0f"]
ARCH_MAP_UNMATCHED_COLOR = "#64748b"


def render_architecture_map_svg(group_rows):
    """group_rows: list of (row_label, motif_list). Renders a SCHEMATIC
    core-promoter-style diagram (same look as a textbook figure):

      * one double-line "track" per group, elements drawn in TSS-relative
        order as equal-size boxes, evenly spaced -> text can never overlap;
      * bp range printed ABOVE each box, the element name INSIDE it, and
        (named core elements only) a wrapped description BELOW it;
      * every element is upstream: the track ends at the red dashed TSS
        line, and an element whose range would pass 0 is clamped to end at 0.

    Because spacing is even, the diagram is NOT to scale -- the printed
    bp ranges carry the true (median, TSS-relative) positions.
    Returns "" if there is nothing to draw.
    """
    from html import escape
    import textwrap

    group_rows = [(label, motifs) for label, motifs in group_rows if motifs]
    if not group_rows:
        return ""

    FONT = "Arial, Helvetica, sans-serif"
    box_h = 46
    gap = 26
    PX_PER_BP = 10            # box width = element length (bp) x this ...
    MIN_BOX_W = 84            # ... but never narrower than this, so names/ranges stay legible
    left_pad = 190
    tss_gap = 46              # track continues this far past the last box, to the TSS
    right_pad = 40
    top_pad = 44
    FS_RANGE, FS_NAME, FS_DESC = 12, 13, 11
    DESC_LINE_H = 14
    ROW_GAP = 26

    def fmt_pos(v):
        return f"+{v}" if v > 0 else str(v)

    label_colors = {}

    def color_for(m):
        if m.get("is_core_element"):
            return color_for_core_element(m["label"])
        if not m["matched"]:
            return ARCH_MAP_UNMATCHED_COLOR
        if m["label"] not in label_colors:
            label_colors[m["label"]] = ARCH_MAP_COLORS[len(label_colors) % len(ARCH_MAP_COLORS)]
        return label_colors[m["label"]]

    def wrap(text, chars, keep_words=False):
        return textwrap.wrap(str(text), chars, break_long_words=not keep_words,
                             break_on_hyphens=not keep_words) or [""]

    # ---- prepare rows: clamp to upstream, sort by position -------------
    rows = []
    for row_label, motifs in group_rows:
        items = []
        for m in motifs:
            # Only annotated elements are drawn: core promoter elements and
            # TOMTOM-matched motifs. Unannotated ("Motif N") motifs are dropped.
            if not m.get("is_core_element") and not m.get("matched"):
                continue
            w = max(int(m["width"]), 1)
            start = int(m["start_tss"])
            end = start + w
            if end > 0:
                start -= end
                end = 0
            desc = ""
            if m.get("is_core_element"):
                d = description_for_core_element(m["label"])
                if d and d.lower() != str(m["label"]).lower():
                    desc = d
            bw = max(MIN_BOX_W, int(round(w * PX_PER_BP)))
            name_lines = wrap(m["label"], 11, keep_words=True)[:3]
            longest = max(len(l) for l in name_lines)
            name_fs = max(9.0, min(float(FS_NAME), (bw - 20) / (longest * 0.6)))
            items.append({"m": m, "start": start, "end": end, "bw": bw,
                          "name_lines": name_lines, "name_fs": name_fs,
                          "desc_lines": wrap(desc, 15) if desc else []})
        items.sort(key=lambda i: (i["start"], i["end"]))
        cursor = left_pad
        for it in items:
            it["x"] = cursor
            cursor += it["bw"] + gap
        if items:
            rows.append((row_label, items))
    if not rows:
        return ""

    plot_w = max(sum(i["bw"] for i in items) + (len(items) - 1) * gap
                 for _, items in rows) + tss_gap
    width = left_pad + plot_w + right_pad

    # ---- row heights ---------------------------------------------------
    row_geo, row_axis, y = [], [], top_pad
    AXIS_GAP = 20            # descriptions -> range baseline
    AXIS_LABEL_H = 22        # tick + range text under the baseline
    for _, items in rows:
        max_desc = max((len(i["desc_lines"]) for i in items), default=0)
        above = 22                                   # range label
        below = 16 + max_desc * DESC_LINE_H if max_desc else 0
        baseline_y = y + above + box_h / 2
        axis_y_row = baseline_y + box_h / 2 + below + AXIS_GAP
        row_geo.append(baseline_y)
        row_axis.append(axis_y_row)
        y = axis_y_row + AXIS_LABEL_H + ROW_GAP
    # ---- legend: which colour = what (wraps onto as many lines as needed) ----
    FS_LEG = 12
    legend_items, seen = [], set()
    for _, items in rows:
        for it in items:
            m = it["m"]
            key = (m["label"], bool(m.get("is_core_element")), bool(m["matched"]))
            if key in seen:
                continue
            seen.add(key)
            if m.get("is_core_element"):
                d = description_for_core_element(m["label"])
                text = f'{m["label"]} ({d})' if d and d.lower() not in str(m["label"]).lower() else str(m["label"])
                legend_items.append((color_for_core_element(m["label"]), text, 1))
            elif m["matched"]:
                legend_items.append((color_for(m), str(m["label"]), 2))
    if any((not it["m"]["matched"] and not it["m"].get("is_core_element")) for _, items in rows for it in items):
        legend_items.append((ARCH_MAP_UNMATCHED_COLOR, "Unmatched motif (no significant TOMTOM hit)", 3))
    legend_items.sort(key=lambda t: t[2])   # core elements, then TF matches, then unmatched

    legend_svg, lx, ly = [], left_pad, y - ROW_GAP + 30
    for color, text, _ in legend_items:
        item_w = 20 + len(text) * FS_LEG * 0.56 + 28
        if lx + item_w > width - 10 and lx > left_pad:
            lx, ly = left_pad, ly + 22
        legend_svg.append(f'<rect x="{lx:.1f}" y="{ly - 11}" width="14" height="14" rx="2" fill="{color}"/>')
        legend_svg.append(f'<text x="{lx + 20:.1f}" y="{ly}" fill="#1e293b" font-size="{FS_LEG}">{escape(text)}</text>')
        lx += item_w
    note_y = (ly + 32) if legend_items else (y - ROW_GAP + 26)
    height = note_y + 12

    tss_x = left_pad + plot_w
    svg = [f'<svg viewBox="0 0 {width} {height:.0f}" xmlns="http://www.w3.org/2000/svg" '
           f'class="arch-map-svg" role="img" font-family="{FONT}">']

    for (row_label, items), baseline_y in zip(rows, row_geo):
        # row label (wrapped, vertically centered on the track)
        rl = wrap(str(row_label), 22)
        ty = baseline_y + 5 - (len(rl) - 1) * 8
        for k, ln in enumerate(rl):
            svg.append(f'<text x="8" y="{ty + k * 16:.1f}" fill="#1e293b" font-size="14">{escape(ln)}</text>')

        # double-line track, ending at the TSS
        for off in (-7, 7):
            svg.append(f'<line x1="{left_pad - 20}" y1="{baseline_y + off:.1f}" x2="{tss_x:.1f}" '
                       f'y2="{baseline_y + off:.1f}" stroke="#0f172a" stroke-width="1.4"/>')

        box_top = baseline_y - box_h / 2
        for idx, it in enumerate(items):
            m = it["m"]
            color = color_for(m)
            x, bw = it["x"], it["bw"]
            cx = x + bw / 2
            pooled_note = " [genus-pooled]" if m.get("genus_pooled") else ""
            svg.append(
                f'<rect x="{x}" y="{box_top:.1f}" width="{bw}" height="{box_h}" rx="{box_h / 2:g}" ry="{box_h / 2:g}" fill="{color}">'
                f'<title>{escape(str(m["label"]))}{pooled_note}'
                f'{(" (MEME " + m["meme_motif"] + ", TOMTOM q=" + format(m["qvalue"], ".2g") + ")") if m.get("matched") and m.get("meme_motif") else ""}'
                f' — median {m["start_tss"]}bp from TSS ({m["n_sites"]} occurrence(s))</title></rect>')
            # range above
            svg.append(f'<text x="{cx:.1f}" y="{box_top - 8:.1f}" fill="#1e293b" font-size="{FS_RANGE}" '
                       f'text-anchor="middle">{fmt_pos(it["start"])} to {fmt_pos(it["end"])}</text>')
            # name inside (1-3 centered lines)
            txt = _contrast_text_color(color)
            nl = it["name_lines"]
            ny = baseline_y + 5 - (len(nl) - 1) * 8
            for k, ln in enumerate(nl):
                svg.append(f'<text x="{cx:.1f}" y="{ny + k * 16:.1f}" fill="{txt}" font-size="{it["name_fs"]:.1f}" '
                           f'font-weight="bold" text-anchor="middle">{escape(ln)}</text>')
            # description below
            for k, ln in enumerate(it["desc_lines"]):
                svg.append(f'<text x="{cx:.1f}" y="{box_top + box_h + 18 + k * DESC_LINE_H:.1f}" '
                           f'fill="#475569" font-size="{FS_DESC}" text-anchor="middle">{escape(ln)}</text>')

    # Range baseline: one per row, directly under the boxes/descriptions. A
    # tick sits under every element's centre with its TSS-relative range
    # printed below it, and the baseline ends at the TSS (0).
    for (row_label, items), axis_y_row in zip(rows, row_axis):
        svg.append(f'<line x1="{left_pad - 20}" y1="{axis_y_row:.1f}" x2="{tss_x:.1f}" y2="{axis_y_row:.1f}" '
                   f'stroke="#0f172a" stroke-width="1.4"/>')
        for idx, it in enumerate(items):
            cx = it["x"] + it["bw"] / 2
            svg.append(f'<line x1="{cx:.1f}" y1="{axis_y_row - 5:.1f}" x2="{cx:.1f}" y2="{axis_y_row:.1f}" '
                       f'stroke="#0f172a" stroke-width="1.4"/>')
            svg.append(f'<text x="{cx:.1f}" y="{axis_y_row + 17:.1f}" fill="#334155" font-size="{FS_RANGE - 1}" '
                       f'text-anchor="middle">{fmt_pos(it["start"])} to {fmt_pos(it["end"])}</text>')
        svg.append(f'<line x1="{tss_x:.1f}" y1="{axis_y_row - 5:.1f}" x2="{tss_x:.1f}" y2="{axis_y_row:.1f}" '
                   f'stroke="#0f172a" stroke-width="1.4"/>')
        svg.append(f'<text x="{tss_x:.1f}" y="{axis_y_row + 17:.1f}" fill="#334155" font-size="{FS_RANGE - 1}" '
                   f'font-weight="bold" text-anchor="middle" stroke="#ffffff" stroke-width="4" paint-order="stroke">0</text>')

    # TSS marker: dashed red line through all rows, label on top
    svg.append(f'<line x1="{tss_x:.1f}" y1="{top_pad - 8}" x2="{tss_x:.1f}" y2="{row_axis[-1]:.1f}" '
               f'stroke="#dc2626" stroke-width="1.6" stroke-dasharray="5,4"/>')
    svg.append(f'<text x="{tss_x:.1f}" y="{top_pad - 16}" fill="#1e293b" font-size="13" '
               f'font-weight="bold" text-anchor="middle">TSS</text>')

    svg.extend(legend_svg)
    svg.append(f'<text x="{left_pad}" y="{note_y:.0f}" fill="#64748b" font-size="10" font-style="italic">'
               f'Upstream of TSS only. Box width is proportional to element length (bp); elements are shown in '
               f'positional order at equal spacing (distances not to scale); labels give TSS-relative bp ranges.</text>')
    svg.append("</svg>")
    return "".join(svg)


def wrap_arch_map_figure(svg, filename_stem):
    """Wraps a rendered inline-SVG architecture map so clicking it opens
    the same zoomable lightbox (with Download + Close) as every other
    panel image in this dashboard -- no separate button, just click the
    map itself, exactly like every base64 <img> panel already works.

    WHY THIS NEEDS ITS OWN HELPER AT ALL: every other panel is a
    base64-encoded <img>, so the lightbox's own download link handles it
    for free from that src. The architecture map is generated as inline
    SVG markup instead (see render_architecture_map_svg), which has no
    src to hand to that link -- so the onclick here rasterizes the live
    SVG on the fly first, then opens it in the normal lightbox.

    The onclick passes `this`, and the JS walks up to `.arch-map-figure`
    to find its own SVG -- so no unique element IDs are needed and
    multiple maps on one page (the species view renders one per genus)
    can't collide."""
    safe_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", filename_stem).strip("_") or "architecture_map"
    return f"""
                <div class="arch-map-figure" onclick="openArchMapLightbox(this, '{safe_stem}.png')">
                    <div class="arch-map-scroll">{svg}</div>
                </div>"""


def render_promoter_architecture_map_section(meme_dir, tomtom_dir, groups, group_by_genus=False):
    """groups: list of (id, label) tuples (genus names, or species ids for
    the species-wise dashboard). group_by_genus bucket-splits the species
    view into one small map per genus with a sub-heading, matching how
    every other panel in the species dashboard already sub-groups.

    Each row combines TWO independent data sources into one set of bars:
      - MEME-discovered motifs matched against TOMTOM (summarize_group_motifs)
      - Named core promoter elements -- BRE/TATA/Inr/etc. -- from
        Core_Promoter_Results_By_Organism's *_individual_hits.csv
        (summarize_core_promoter_elements_for_group), which already gives
        exact TSS-relative positions directly, no inference needed.
    is_genus is inferred from which meme_dir was passed, since that's
    already how the caller distinguishes genus-wise vs species-wise here.
    """
    is_genus = (meme_dir == MEME_OUTPUT_DIR_GENUS)
    normalized = [g if isinstance(g, tuple) else (g, g) for g in groups]

    # Tracks whether ANY rendered bar came from the genus-pooled fallback,
    # so the caption can say so. Set during the single build pass below --
    # build_row() is deliberately called exactly once per group (it does
    # real file I/O), rather than re-running it later just to inspect flags.
    pooled_used = {"value": False}

    def build_row(gid, label):
        bars = summarize_group_motifs(meme_dir, tomtom_dir, gid)
        bars += summarize_core_promoter_elements_for_group(gid, is_genus)
        if any(m.get("genus_pooled") for m in bars):
            pooled_used["value"] = True
        return label, bars

    if group_by_genus:
        genus_buckets = {}
        for gid, label in normalized:
            genus = str(gid).split("_")[0]
            genus_buckets.setdefault(genus, []).append((gid, label))

        sub_html = []
        any_data = False
        for genus, members in genus_buckets.items():
            rows = [build_row(gid, label) for gid, label in members]
            svg = render_architecture_map_svg(rows)
            if not svg:
                continue
            any_data = True
            sub_html.append(f"""
                <div class="genus-subgroup">
                    <div class="genus-subheading">{genus}</div>
                    {wrap_arch_map_figure(svg, f"promoter_architecture_map_{genus}")}
                </div>""")
        if not any_data:
            return ""
        body = "".join(sub_html)
    else:
        rows = [build_row(gid, label) for gid, label in normalized]
        svg = render_architecture_map_svg(rows)
        if not svg:
            return ""
        body = wrap_arch_map_figure(svg, "promoter_architecture_map")

    return f"""
        <div class="section arch-map-section">
            <h2 class="section-title">Promoter Architecture Map</h2>
            {body}
        </div>"""


def _render_one_tree_card(img_path, card_title=None):
    """Builds one <img> card (with staleness banner) for a single tree
    image path. Shared by the combined-tree and per-genus-tree cases below
    so both get identical staleness-checking and lightbox behavior."""
    b64_str = image_to_base64(img_path)
    if not b64_str:
        return ""

    filename = os.path.basename(img_path)
    js_safe_filename = filename.replace("\\", "\\\\").replace("'", "\\'")

    stale_warning_html = ""
    try:
        tree_mtime = os.path.getmtime(img_path)
        labels_mtime = os.path.getmtime(LABELS_PATH) if os.path.exists(LABELS_PATH) else None
        if labels_mtime is not None and labels_mtime > tree_mtime:
            stale_warning_html = """
            <div class="stale-warning" style="background:#4a2a0a; border:1px solid #d98c2b;
                 color:#ffd9a0; padding:10px 14px; border-radius:6px; margin-bottom:10px;
                 font-size:0.9em;">
                &#9888; <strong>This tree may be out of date.</strong>
                family_labels.csv has been modified more recently than this image, which is
                built once by the OrthoFinder/STAG step (03c) and not regenerated
                automatically. If you've changed the organism set since the tree was last
                built, re-run 03c before trusting this panel.
            </div>"""
    except OSError:
        pass  # if mtimes aren't readable for any reason, skip the check rather than block rendering

    title_html = f'<div class="reference-tree-card-title">{card_title}</div>' if card_title else ""

    return f"""
            <div class="reference-tree-card">
                {title_html}
                {stale_warning_html}
                <img src="{b64_str}" alt="Reference species tree"
                     onclick="openLightbox('{b64_str}', '{js_safe_filename}')"/>
            </div>"""


def render_reference_species_tree_section(genera=None):
    """Single, non-repeated section shown once at the top of BOTH
    dashboards -- the OrthoFinder/STAG reference species tree isn't a
    per-genus or per-species comparison like every other panel, it's one
    tree answering "how do the host organisms relate", so it doesn't fit
    the group_word-vs-group_word card-row pattern the rest of this file
    uses. Returns "" (renders nothing) if no tree has been built yet at
    all -- e.g. OrthoFinder isn't installed, or 03c's reference-tree step
    hasn't run -- rather than showing a permanent "Missing" card for
    something that may be an optional step for this run.

    03c attempts TWO kinds of tree: one "all genera combined" tree, and
    (separately) one tree PER GENUS. Either can succeed or fail
    independently -- in practice, a genus with too few resulting
    orthogroups (e.g. Babesia with 3 very divergent species) can fail its
    species tree while another genus (e.g. Toxoplasma) succeeds, and the
    combined-all-genera run can *also* fail even when a per-genus run
    didn't. This function therefore shows whichever combined and/or
    per-genus tree images actually exist on disk, instead of only ever
    looking for the single combined-tree path -- a genus's tree having
    built successfully shouldn't be hidden just because a different run
    (or the combined run) didn't produce one.

    STALENESS CHECK: each tree image found is a static file built once by
    an earlier pipeline stage (03c) and not regenerated automatically when
    the rest of the dashboard's data changes. Every other panel here is
    re-derived from the CURRENT family_labels.csv each time this script
    runs, but these tree images are not -- if family_labels.csv is later
    regenerated with a different organism set without re-running 03c,
    the image still holds the OLD tree. The best available proxy is
    comparing file modification times: if family_labels.csv is newer than
    a given tree image, that tree almost certainly predates the current
    organism set. That's flagged with a visible on-page banner per card
    (not just a console print) since a silently wrong tree is easy to
    miss otherwise -- the image is still shown rather than hidden, since a
    stale tree is still more informative than no tree at all."""
    cards_html = []

    combined_matches = find_images_for_group("", REFERENCE_SPECIES_TREE_PATTERNS, OUTPUT_DIR, False)
    if combined_matches:
        title = "All genera combined" if genera and len(genera) > 1 else None
        card = _render_one_tree_card(combined_matches[0], card_title=title)
        if card:
            cards_html.append(card)

    for genus in (genera or []):
        genus_matches = find_images_for_group(
            genus, REFERENCE_SPECIES_TREE_PATTERNS_BY_GENUS, OUTPUT_DIR, False
        )
        if genus_matches:
            card = _render_one_tree_card(genus_matches[0], card_title=f"Genus: {genus}")
            if card:
                cards_html.append(card)

    if not cards_html:
        return ""

    return f"""
        <div class="section reference-tree-section">
            <h2 class="section-title">Reference Species Tree (OrthoFinder / STAG)</h2>
            <div class="reference-tree-card-row" style="display:flex; flex-wrap:wrap; gap:16px;">
                {''.join(cards_html)}
            </div>
        </div>"""


def render_full_document(page_title, group_word, groups, sections, header_actions_html,
                          group_by_genus=False, meme_dir=None, tomtom_dir=None, genera=None):
    """Wraps rendered sections in a COMPLETE, standalone HTML document (own
    <html>/<head>/CSS/lightbox script). Used both for the genus-wise file
    actually written to disk, and for the species-wise document that gets
    embedded (as a JS string, never as its own file) inside it."""
    normalized_groups = [g if isinstance(g, tuple) else (g, g) for g in groups]
    group_display_list = ", ".join(label for _, label in normalized_groups)
    sections_html, n_found, n_missing = render_sections_html(sections, groups, group_by_genus)
    kpi_html = render_kpi_cards_html()
    reference_tree_html = render_reference_species_tree_section(genera)
    arch_map_html = ""
    if meme_dir and tomtom_dir:
        arch_map_html = render_promoter_architecture_map_section(
            meme_dir, tomtom_dir, groups, group_by_genus
        )

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{page_title}</title>
    <style>{CSS_BLOCK}</style>
</head>
<body>
    <div class="dashboard-container">
        <header>
            <div class="header-text">
                <h1>{page_title}</h1>
                <p class="subtitle">{group_word}-vs-{group_word} comparison — {group_display_list} (click any image to zoom/download)</p>
            </div>
            {header_actions_html}
        </header>
        {kpi_html}
        <div class="sections-stack">{reference_tree_html}{arch_map_html}{sections_html}
        </div>
    </div>
{LIGHTBOX_HTML_AND_JS}
</body>
</html>
"""
    return html_content, n_found, n_missing


def main():
    print("Building comparative dashboard (genus-wise, with species-wise "
          "available as an on-the-fly download)...")

    genera = GENERA_OVERRIDE if GENERA_OVERRIDE else detect_genera(LABELS_PATH)
    species_ids = SPECIES_OVERRIDE if SPECIES_OVERRIDE else detect_species(LABELS_PATH)
    # Nicer display label for species cards: underscores back to spaces.
    species_groups = [(sid, sid.replace("_", " ")) for sid in species_ids]

    if not genera and not species_ids:
        print("[!] No genera or species detected/configured -- nothing to compare. "
              "Set GENERA_OVERRIDE/SPECIES_OVERRIDE explicitly at the top of this script.")
        return

    # Build the species-wise document FIRST. It's never written to its own
    # file -- its rendered HTML is only ever used as the string the download
    # button hands to the browser. header_actions_html="" here: the
    # downloaded copy has no further download button of its own (a
    # reciprocal "download genus-wise" button would mean embedding a SECOND
    # full document inside this one, which isn't worth the complexity for
    # what's meant to stay a light on-the-fly export).
    species_html, species_found, species_missing = "", 0, 0
    if species_ids:
        species_html, species_found, species_missing = render_full_document(
            "Comparative Promoter Pipeline Dashboard — Species-wise",
            "Species", species_groups, SPECIES_SECTIONS,
            header_actions_html="",
            group_by_genus=True,
            meme_dir=MEME_OUTPUT_DIR_ORGANISM,
            tomtom_dir=TOMTOM_OUTPUT_DIR_ORGANISM,
            genera=genera,
        )
        print(f"[*] Species-wise data embedded for on-the-fly download: "
              f"{', '.join(species_ids)} ({species_found} found, {species_missing} missing)")
    else:
        print("[!] No species detected -- the download button will be omitted "
              "from the dashboard.")

    if not genera:
        print("[!] No genera detected -- cannot build the main (genus-wise) dashboard.")
        return

    header_actions_html = DOWNLOAD_BUTTON_HTML if species_html else ""

    genus_html, genus_found, genus_missing = render_full_document(
        "Comparative Promoter Pipeline Dashboard — Genus-wise",
        "Genus", genera, GENUS_SECTIONS,
        header_actions_html=header_actions_html,
        meme_dir=MEME_OUTPUT_DIR_GENUS,
        tomtom_dir=TOMTOM_OUTPUT_DIR_GENUS,
        genera=genera,
    )

    if species_html:
        # See module docstring's "CRITICAL SUBTLETY" section for why the
        # "</script" replace is mandatory, not optional, given the embedded
        # document contains its own <script> block.
        escaped_species_json = json.dumps(species_html).replace("</script", "<\\/script")
        download_script = DOWNLOAD_SCRIPT_TEMPLATE.format(
            species_html_json=escaped_species_json,
            species_filename_json=json.dumps(SPECIES_DOWNLOAD_FILENAME),
        )
        genus_html = genus_html.replace("</body>", f"{download_script}</body>")

    with open(DASHBOARD_FILE, "w", encoding="utf-8") as f:
        f.write(genus_html)

    print(f"[*] Genus-wise dashboard comparing: {', '.join(genera)}")
    print(f"    -> {DASHBOARD_FILE}  ({genus_found} found, {genus_missing} missing)")
    if species_html:
        print(f"    Top-right button downloads the species-wise dashboard on the fly "
              f"as '{SPECIES_DOWNLOAD_FILENAME}' -- no second file is written to disk "
              f"by this script.")

    print(
        "\n[!] Reminder: several panel patterns above are educated guesses, not "
        "confirmed against every stage script's exact output filename. A panel "
        "showing 'Missing' across the board is the signal to check that panel's "
        "'patterns' list against your actual outputs/ tree."
    )


if __name__ == "__main__":
    main()
