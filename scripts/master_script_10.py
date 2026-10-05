#!/usr/bin/env python3

"""
==============================================================================
        COMPARATIVE GENE FAMILY & PROMOTER ANALYSIS
                        MASTER PIPELINE
==============================================================================

Author: Swarup Das
        Subarna Thakur
        Department of Bioinformatics, University of North Bengal
        Bairatisal, West Bengal 734013
Purpose:
    Master controller for the comparative gene-family and promoter-analysis
    pipeline.

IMPORTANT:
    Stage 01 is intentionally INTERACTIVE.

    Stage 01:
        01e_new3_fetch_family_and_promoters.py

    The master script does NOT collect Stage 01 input itself.

    Instead, Stage 01 is launched with inherited stdin/stdout/stderr so that
    its interactive prompts behave exactly as they do when the user runs:

        python 01e_new3_fetch_family_and_promoters.py

    After Stage 01 finishes successfully, all subsequent stages are run
    automatically.

==============================================================================

PIPELINE STRUCTURE

    Stage 01  -> Interactive data retrieval
    Stage 02  -> Phylogenetic tree
    Stage 03  -> MSA / conservation
    Stage 04  -> TSS prediction
    Stage 05  -> Nucleotide composition
    Stage 06  -> Core promoter elements
    Stage 07  -> Core promoter synergism
    Stage 08  -> Motif discovery
    Stage 09  -> Motif logos
    Stage 10  -> Motif conservation
    Stage 11  -> Motif distribution
    Stage 12  -> CpG islands
    Stage 13  -> DNA structural properties
    Stage 14  -> Non-B DNA motifs

    MEME-related stages use the MEME Suite executables available on PATH.

==============================================================================

REQUIREMENTS

    Python
    MEME Suite tools available on PATH

MEME tools checked:
    meme
    tomtom
    ceqlogo

MOTIF REFERENCE DATABASES:
    JASPAR2026_CORE.meme and UniProbe_Combined.meme are set up ONCE into a
    fixed databases/ folder next to this script (via
    00_setup_motif_databases.py, invoked automatically as part of
    pre-flight checks below) rather than being re-downloaded into each
    run's outputs/ folder. Needs network access the first time only.

==============================================================================
"""

import os
import sys
import time
import json
import argparse
import subprocess
import shutil
from pathlib import Path
from datetime import datetime


# =============================================================================
# CONFIGURATION
# =============================================================================

# SCRIPT_DIR : where the stage scripts live. When the pipeline is installed as a conda
#              package this is inside the environment and must be treated as read-only.
# WORK_DIR   : where the user's input data is and where outputs/ is written. It is the
#              directory the pipeline is started from (or COMPARATIVE_PIPELINE_WORKDIR).
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR  # backward-compatible alias; stage scripts are found here
WORK_DIR = Path(os.environ.get("COMPARATIVE_PIPELINE_WORKDIR") or Path.cwd()).resolve()

_conda_prefix = os.environ.get("CONDA_PREFIX")
if _conda_prefix:
    try:
        WORK_DIR.relative_to(Path(_conda_prefix).resolve())
    except ValueError:
        pass  # normal case: working folder is outside the conda environment
    else:
        sys.exit(
            "[ERROR] The working directory is inside the conda environment:\n"
            f"        {WORK_DIR}\n"
            "        Outputs would be written into the installed package. "
            "cd to your own data folder first."
        )

OUTPUT_DIR = WORK_DIR / "outputs"
LOG_DIR = OUTPUT_DIR / "logs"
STATE_FILE = OUTPUT_DIR / ".pipeline_state.json"
DATA_MANIFEST_PATH = OUTPUT_DIR / ".pipeline_data_manifest.json"
# When 01e_new4_fetch_family_and_promoters.py detects no real protein data
# for this run (local genome+GFF mode extracts nucleotide gene sequences,
# not translated protein -- see that script's detect_sequence_type()), this
# pipeline enforces starting at Stage NO_PROTEIN_MIN_STAGE and skips any
# stage marked "requires_protein": True in STAGES, wherever it sits.
NO_PROTEIN_MIN_STAGE = 5


# -----------------------------------------------------------------------------
# Python executable
# -----------------------------------------------------------------------------

PYTHON = sys.executable


# =============================================================================
# STAGE DEFINITIONS
# =============================================================================

STAGES = [

    # -------------------------------------------------------------------------
    # STAGE 01
    # IMPORTANT:
    # This stage MUST remain interactive.
    # -------------------------------------------------------------------------

    {
        "number": 1,
        "name": "Gene Family, Domain & Promoter Data Extraction",
        "script": "01e_new9_fetch_family_and_promoters.py",
        "environment": "current",
        "interactive": True,
    },


    # -------------------------------------------------------------------------
    # STAGE 02
    # Only needs outputs/family_sequences.fasta (+ outputs/family_labels.csv
    # for the organism/genus breakdown) from Stage 01 -- no dependency on
    # alignment or phylogeny, so it can run immediately after retrieval.
    # 02b imports 02a_protein_properties_for_all.py directly (same
    # base-script/wrapper convention as 05f importing 05c/05d) and also
    # reproduces its global CSV/plot, so only 02b needs to be listed here.
    # -------------------------------------------------------------------------

    {
        "number": 2,
        "name": "Protein Physicochemical Property Analysis (organism/genus-wise)",
        "script": "02b_organism_wise_protein_properties_export.py",
        "environment": "current",
        "interactive": False,
        "requires_protein": True,
    },


    # -------------------------------------------------------------------------
    # STAGE 03
    # Depends on outputs/family_aligned.fasta produced by Stage 02 above.
    # -------------------------------------------------------------------------

    {
        "number": 3,
        "name": "Neighbor-Joining Phylogenetic Tree ",
        "script": "03c_new_seperate_orthofinder_run_organism_wise_phylogenetic_trees.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 04
    # Also depends on outputs/family_aligned.fasta produced by Stage 02 above.
    # -------------------------------------------------------------------------

    {
        "number": 4,
        "name": "Multiple Sequence Alignment & Conservation Profiling",
        "script": "04b_new_new_new_organism_wise_msa_visualization.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 05
    # -------------------------------------------------------------------------

    {
        "number": 5,
        "name": "TSS Location Prediction",
        "script": "05g_tss_prediction_upstream_only_single_script.py",
        "environment": "current",
        "interactive": False,
    },


    # -------------------------------------------------------------------------
    # STAGE 06
    # -------------------------------------------------------------------------

    {
        "number": 6,
        "name": "Nucleotide Frequency & Windowed Composition",
        "script": "06a_new_organism_wise_nucleotide_composition_analysis.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 08
    # -------------------------------------------------------------------------

    {
        "number": 7,
        "name": "Core Promoter Element Analysis",
        "script": "06b_new_organism_wise_core_promoter_analysis.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 09
    # -------------------------------------------------------------------------

    {
        "number": 8,
        "name": "Core Promoter Element Synergism / Co-occurrence",
        "script": "06c_new_organism_wise_synergism_analysis.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 10
    # MEME ENVIRONMENT
    # -------------------------------------------------------------------------

    {
        "number": 9,
        "name": "Motif Discovery",
        "script": "07c_new1_Motif_discovery_for_all.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 11
    # MEME ENVIRONMENT
    # -------------------------------------------------------------------------

    {
        "number": 10,
        "name": "Motif Logo Generation",
        "script": "07e_new_generate_motif_logos_for_all.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 12
    # MEME ENVIRONMENT
    # -------------------------------------------------------------------------

    {
        "number": 11,
        "name": "Motif Conservation Analysis",
        "script": "07f_new5_for_all_motif_conservation_analysis.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 13
    # -------------------------------------------------------------------------

    {
        "number": 12,
        "name": "Motif Distribution Analysis",
        "script": "07i_new_organism_and_genus_motif_distribution.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 14
    # -------------------------------------------------------------------------

    {
        "number": 13,
        "name": "CpG Island Analysis",
        "script": "08_new_organism_wise_cpg_island_finder.py",
        "environment": "current",
        "interactive": False,
    },

# -------------------------------------------------------------------------
    # STAGE 14
    # -------------------------------------------------------------------------

    {
        "number": 14,
        "name": "CpG Island Analysis Genewise",
        "script": "08c_gene_wise_cpg_methprimer_plots.py",
        "environment": "current",
        "interactive": False,
    },
    
    # -------------------------------------------------------------------------
    # STAGE 15
    # -------------------------------------------------------------------------

    {
        "number": 15,
        "name": "DNA Structural Properties",
        "script": "09_new4_dna_structural_properties_pipeline.py",
        "environment": "current",
        "interactive": False,
    },

    # -------------------------------------------------------------------------
    # STAGE 16
    # -------------------------------------------------------------------------

    {
        "number": 16,
        "name": "Non-B DNA Motif Analysis",
        "script": "010_new1_nonB_dna_motif_analysis_for_all.py",
        "environment": "current",
        "interactive": False,
    },
    
     # -------------------------------------------------------------------------
    # STAGE 17
    # -------------------------------------------------------------------------

    {
        "number": 17,
        "name": "Dashboard generation",
        "script": "011u5u_generate_comparative_dashboard.py",
        "environment": "current",
        "interactive": False,
    },
]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def timestamp():
    """Return current timestamp."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def print_header():
    """Print main pipeline header."""

    print()
    print("=" * 78)
    print("        COMPARATIVE GENE FAMILY & PROMOTER ANALYSIS")
    print("=" * 78)
    print()
    print("  Developed by: Swarup Das & Subarna Thakur")
    print("  Department of Bioinformatics, University of North Bengal")
    print("  Raja Rammohunpur, Bagdogra, Bairatisal, West Bengal 734013")
    print("Working directory:", WORK_DIR)
    print("Output directory :", OUTPUT_DIR)
    print("Python executable:", PYTHON)
    print("Conda environment:", os.environ.get("CONDA_DEFAULT_ENV", "not activated"))
    print("MEME tools       : PATH")
    print()
    print("Started:", timestamp())
    print("=" * 78)
    print()


def ensure_directories():
    """Create required output and log directories."""

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    LOG_DIR.mkdir(
        parents=True,
        exist_ok=True
    )


# =============================================================================
# RESUME / PAUSE STATE
#
# Every time a stage finishes successfully, its number is appended to
# STATE_FILE (outputs/.pipeline_state.json). This lets a later run pick up
# with --resume (continue after the last successful stage) or --from N
# (jump straight to stage N, e.g. "after retrieval" or "after phylogeny"),
# without re-running everything from scratch.
# =============================================================================

def load_state():
    """
    Read which stages have completed successfully in previous runs.
    Never raises -- a missing or corrupt state file just means "nothing
    completed yet" rather than crashing the pipeline.
    """

    if not STATE_FILE.exists():
        return {"completed_stages": []}

    try:
        with open(STATE_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
        if "completed_stages" not in data:
            data["completed_stages"] = []
        return data
    except Exception as e:
        print(f"[WARN] Could not read state file ({e}) -- treating as no progress saved.")
        return {"completed_stages": []}


def save_state(stage_number):
    """Record that a stage finished successfully, for future --resume runs."""

    ensure_directories()

    data = load_state()
    completed = set(data.get("completed_stages", []))
    completed.add(stage_number)

    data["completed_stages"] = sorted(completed)
    data["last_updated"] = timestamp()

    try:
        with open(STATE_FILE, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    except Exception as e:
        # Never let a state-file write failure take down an otherwise
        # successful pipeline run -- just warn about it.
        print(
            f"[WARN] Could not write state file ({e}). "
            f"--resume may not reflect this stage next time."
        )


def reset_state():
    """Delete saved progress so the next run starts completely fresh."""

    if STATE_FILE.exists():
        STATE_FILE.unlink()
        print(f"[OK] Cleared saved progress ({STATE_FILE}).")
    else:
        print("[OK] No saved progress to clear.")


def get_last_completed_stage():
    """Highest stage number recorded as successfully completed, or None."""

    completed = load_state().get("completed_stages", [])
    return max(completed) if completed else None


def load_data_manifest():
    """
    Read outputs/.pipeline_data_manifest.json, written by
    01e_new4_fetch_family_and_promoters.py's write_data_manifest() after
    Stage 01. Returns a dict; if the manifest is missing (e.g. an older
    fetcher script that doesn't write one, or Stage 01 hasn't run yet in
    this outputs/ directory), has_real_protein defaults to True so runs
    that never touch protein data are completely unaffected.
    """

    if not DATA_MANIFEST_PATH.exists():
        return {"has_real_protein": True}

    try:
        with open(DATA_MANIFEST_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        if "has_real_protein" not in data:
            data["has_real_protein"] = True
        return data
    except Exception as e:
        print(
            f"[WARN] Could not read data manifest ({e}) -- assuming real "
            f"protein data IS available (safest default: nothing gets "
            f"silently skipped)."
        )
        return {"has_real_protein": True}


def resolve_start_stage(requested_number):
    """
    Return the STAGES list index of the first stage whose number is
    >= requested_number.

    Stage numbers in this pipeline aren't always contiguous (a stage can be
    merged/removed, leaving a gap), so "--from 2" on a pipeline with no
    Stage 02 maps onto the next stage that actually exists instead of
    failing outright.
    """

    available = sorted(s["number"] for s in STAGES)

    if requested_number in available:
        target = requested_number
    else:
        later = [n for n in available if n > requested_number]
        if not later:
            print(
                f"[ERROR] No stage with number >= {requested_number} exists. "
                f"Highest stage is {available[-1]:02d}."
            )
            sys.exit(1)
        target = later[0]
        print(
            f"[NOTE] Stage {requested_number:02d} doesn't exist in this "
            f"pipeline -- starting from the next available stage, "
            f"{target:02d}, instead."
        )

    for i, stage in enumerate(STAGES):
        if stage["number"] == target:
            return i

    sys.exit(1)  # unreachable given the checks above


def print_stage_list():
    """Print every stage with its completion status, then exit."""

    completed = set(load_state().get("completed_stages", []))

    print()
    print("=" * 78)
    print("PIPELINE STAGES")
    print("=" * 78)
    print()

    for stage in STAGES:
        mark = "[DONE]" if stage["number"] in completed else "[    ]"
        interactive = "  (interactive)" if stage.get("interactive", False) else ""
        print(f"{mark} Stage {stage['number']:02d}: {stage['name']}{interactive}")
        print(f"           {stage['script']}")

    print()
    print(f"State file: {STATE_FILE}")
    print()
    print("Resume after the last completed stage:")
    print("    python master_pipeline.py --resume")
    print()
    print("Jump straight to a specific stage (e.g. after phylogeny):")
    print("    python master_pipeline.py --from <stage_number>")
    print("=" * 78)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Comparative gene-family and promoter-analysis pipeline."
    )

    start_group = parser.add_mutually_exclusive_group()

    start_group.add_argument(
        "--resume", action="store_true",
        help="Resume automatically after the last successfully completed stage.",
    )
    start_group.add_argument(
        "--from", dest="from_stage", type=int, metavar="N",
        help="Start from stage N, skipping everything before it "
             "(Stage 01's interactive retrieval is skipped and assumed done "
             "whenever N > 1).",
    )

    parser.add_argument(
        "--list", action="store_true",
        help="List all stages with completion status, then exit.",
    )
    parser.add_argument(
        "--reset", action="store_true",
        help="Clear saved resume progress, then exit.",
    )

    return parser.parse_args()


def print_stage_header(stage):
    """Print stage information."""

    print()
    print("=" * 78)
    print(
        f"STAGE {stage['number']:02d}: "
        f"{stage['name']}"
    )
    print("=" * 78)

    print()
    print("Script      :", stage["script"])
    print("Environment :", stage["environment"])

    if stage.get("interactive", False):
        print("Mode        : INTERACTIVE")
    else:
        print("Mode        : AUTOMATIC")

    print()


def get_log_file(stage):
    """Return stage log filename and ensure parent log directory exists."""

    LOG_DIR.mkdir(parents=True, exist_ok=True)

    safe_name = (
        stage["script"]
        .replace(".py", "")
        .replace("/", "_")
        .replace(" ", "_")
    )

    return LOG_DIR / f"stage_{stage['number']:02d}_{safe_name}.log"


# =============================================================================
# ENVIRONMENT CHECKS
# =============================================================================

def check_runtime_environment():
    """Check the active Python interpreter and required external tools.

    The packaged pipeline intentionally does not depend on a separate Conda
    environment name. The executable that launches this master script is the
    Python interpreter used for every Python stage, while MEME/OrthoFinder and
    their command-line dependencies are resolved directly from PATH.
    """

    print("=" * 78)
    print("CHECKING RUNTIME ENVIRONMENT")
    print("=" * 78)
    print()

    print(f"[OK] Python: {PYTHON}")
    print(f"     Version: {sys.version.split()[0]}")
    print(
        f"     Conda environment: "
        f"{os.environ.get('CONDA_DEFAULT_ENV', 'not activated')}"
    )

    tools = [
        "datasets",
        "meme",
        "tomtom",
        "ceqlogo",
        "mafft",
        "orthofinder",
        "diamond",
        "mcl",
        "FastTree",
    ]

    all_ok = True

    print()
    print("Checking required command-line tools:")

    aliases = {"FastTree": ("FastTree", "fasttree", "FastTreeMP")}

    for tool in tools:
        path = next(
            (p for p in (shutil.which(n) for n in aliases.get(tool, (tool,))) if p),
            None,
        )

        if path:
            print(f"[OK] {tool}: {path}")
        else:
            print(f"[ERROR] {tool} not found on PATH")
            all_ok = False

    return all_ok

# =============================================================================
# MEME TOOL CHECK
# =============================================================================

def check_meme_tools():
    """Verify required MEME Suite commands are directly available on PATH."""

    print()
    print("=" * 78)
    print("CHECKING MEME SUITE")
    print("=" * 78)

    tools = ["meme", "tomtom", "ceqlogo"]
    all_ok = True

    for tool in tools:
        path = shutil.which(tool)
        if path:
            print(f"[OK] {tool}")
            print(f"     {path}")
        else:
            print(f"[ERROR] {tool} not found on PATH.")
            all_ok = False

    return all_ok



# =============================================================================
# SCRIPT CHECK
# =============================================================================

def check_pipeline_scripts():
    """
    Check whether all pipeline scripts exist.
    """

    print()
    print("=" * 78)
    print("CHECKING PIPELINE SCRIPTS")
    print("=" * 78)

    missing = []

    for stage in STAGES:

        script_path = SCRIPT_DIR / stage["script"]

        if script_path.exists():

            print(
                f"[OK] Stage {stage['number']:02d}: "
                f"{stage['script']}"
            )

        else:

            print(
                f"[MISSING] Stage {stage['number']:02d}: "
                f"{stage['script']}"
            )

            missing.append(stage["script"])

    print()

    if missing:

        print("[ERROR] Missing pipeline scripts:")

        for script in missing:
            print("   -", script)

        return False

    print("[OK] All primary pipeline scripts found.")

    return True


# =============================================================================
# UNDERLYING SCRIPT CHECK
# =============================================================================

def check_underlying_scripts():
    """
    Check scripts that are used by wrapper/analysis stages.
    """

    print()
    print("=" * 78)
    print("CHECKING UNDERLYING PIPELINE SCRIPTS")
    print("=" * 78)

    underlying_scripts = [

        "02a_protein_properties_for_all.py",

        "02_align_with_mafft_for_all.py",

        "03_Build_neighbor_joining_phylogenetic_tree_with_thousand_bootstraps.py",

        "04_Multiple_sequence_alignment_visualization_conservation_profiling.py",

        "06a_Nucleotide_Frequency_Windowed_Composition_Analysis.py",

        "06b_correct_with_adjustments_for_all_Core_Promoter_Element_Analysis.py",

        "06c_Core_Promoter_Element_Synergism_Co_occurrence_Analysis_for_all.py",

        "08_CpG_Island_adjusted_Finder_for_all.py",

        "07i_motif_distribution_block_diagram.py",

    ]

    missing = []

    for script in underlying_scripts:

        path = SCRIPT_DIR / script

        if path.exists():

            print(f"[OK] {script}")

        else:

            print(f"[MISSING] {script}")

            missing.append(script)

    print()

    if missing:

        print("[ERROR] Underlying scripts are missing.")

        return False

    print("[OK] dependencies found.")

    return True


# =============================================================================
# MOTIF REFERENCE DATABASE CHECK
# =============================================================================

def check_motif_databases():
    """Run the packaged motif-database setup script with the active Python.

    The setup script is responsible for placing JASPAR/UniProbe files in the
    package-level ``databases/`` directory.  It is intentionally invoked on
    every startup because it is idempotent.
    """

    print()
    print("=" * 78)
    print("CHECKING MOTIF REFERENCE DATABASES")
    print("=" * 78)

    setup_script = SCRIPT_DIR / "00_setup_motif_databases.py"

    if not setup_script.exists():
        print(f"[ERROR] {setup_script.name} not found in {SCRIPT_DIR}.")
        return False

    command = [PYTHON, "-u", str(setup_script)]

    try:
        result = subprocess.run(command, cwd=str(WORK_DIR))
    except Exception as e:
        print(f"[ERROR] Could not run {setup_script.name}: {e}")
        return False

    if result.returncode != 0:
        print()
        print(
            f"[ERROR] {setup_script.name} exited with code {result.returncode} -- "
            "motif reference databases are not ready. Motif Discovery and "
            "downstream motif analyses will fail until this is fixed "
            "(network access may be required on the first run)."
        )
        return False

    print()
    print("[OK] Motif reference databases ready.")
    return True

def check_orthofinder_tools():
    """Verify OrthoFinder and its command-line dependencies on PATH.

    This remains informational/non-blocking because the reference species-tree
    enhancement in Stage 03 can skip cleanly when OrthoFinder is unavailable.
    """

    print()
    print("=" * 78)
    print("CHECKING ORTHOFINDER (used by the reference species tree step in Stage 3)")
    print("=" * 78)

    tools = ["orthofinder", "diamond", "mcl"]
    all_ok = True

    for tool in tools:
        path = shutil.which(tool)
        if path:
            print(f"[OK] {tool}")
            print(f"     {path}")
        else:
            print(f"[WARN] {tool} not found on PATH.")
            all_ok = False

    if not all_ok:
        print()
        print("[NOTE] The reference species tree step in Stage 3 will skip itself "
              "cleanly rather than fail; other stages are unaffected.")

    return all_ok


# =============================================================================
# STAGE 01 -- INTERACTIVE
# =============================================================================

def run_stage_01_interactive(stage):
    """
    Run Stage 01 interactively with direct terminal I/O streaming and unbuffered output.
    """

    script_path = SCRIPT_DIR / stage["script"]

    print()
    print("=" * 78)
    print("STAGE 01 IS INTERACTIVE")
    print("=" * 78)

    print()
    print("Stage 01 will now run directly in your terminal.")
    print()
    print("Please answer the prompts displayed by:")
    print(f"    {stage['script']}")
    print()
    print("The master script will NOT collect or modify your Stage 01 inputs.")
    print()
    print("=" * 78)
    print("STARTING STAGE 01")
    print("=" * 78)
    print()

    try:
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"

        # Direct execution is intentional: inherited stdin/stdout/stderr keep
        # Stage 01 fully interactive and avoid conda-run stdin issues.
        command = [PYTHON, "-u", str(script_path)]

        result = subprocess.run(
            command,
            cwd=str(WORK_DIR),
            stdin=sys.stdin,
            stdout=sys.stdout,
            stderr=sys.stderr,
            env=env
        )

        print()

        if result.returncode == 0:
            print("=" * 78)
            print("STAGE 01 COMPLETED SUCCESSFULLY")
            print("=" * 78)
            return True
        else:
            print("=" * 78)
            print(f"STAGE 01 FAILED (exit code: {result.returncode})")
            print("=" * 78)
            return False

    except KeyboardInterrupt:
        print()
        print()
        print("=" * 78)
        print("STAGE 01 INTERRUPTED BY USER")
        print("=" * 78)
        return False

    except Exception as e:
        print()
        print(f"[ERROR] Failed to execute Stage 01: {e}")
        return False


# =============================================================================
# AUTOMATIC STAGES
# =============================================================================

def run_automatic_stage(stage):
    """
    Run all non-interactive stages.

    Output is simultaneously displayed and written to a log file.
    """

    script_path = SCRIPT_DIR / stage["script"]

    log_file = get_log_file(stage)

    print_stage_header(stage)

    print(
        f"Log file: {log_file}"
    )

    print()

    # -------------------------------------------------------------------------
    # Determine command
    # -------------------------------------------------------------------------

    # All Python stages use the same interpreter that launched this master
    # script. MEME tools are resolved from PATH by the stage scripts.
    command = [PYTHON, "-u", str(script_path)]

    # -------------------------------------------------------------------------
    # Start process
    # -------------------------------------------------------------------------

    try:

        with open(
            log_file,
            "w",
            encoding="utf-8"
        ) as log:

            log.write(
                "=" * 78 + "\n"
            )

            log.write(
                f"STAGE {stage['number']:02d}\n"
            )

            log.write(
                f"Script: {stage['script']}\n"
            )

            log.write(
                f"Started: {timestamp()}\n"
            )

            log.write(
                "=" * 78 + "\n\n"
            )

            log.flush()

            process = subprocess.Popen(

                command,

                cwd=str(WORK_DIR),

                stdin=subprocess.DEVNULL,

                stdout=subprocess.PIPE,

                stderr=subprocess.STDOUT,

                text=True,

                bufsize=1
            )

            # -------------------------------------------------------------
            # Stream output to terminal and log simultaneously
            # -------------------------------------------------------------

            if process.stdout is not None:

                for line in process.stdout:

                    print(
                        line,
                        end=""
                    )

                    log.write(line)

                    log.flush()

            return_code = process.wait()

            # -------------------------------------------------------------
            # Write completion information
            # -------------------------------------------------------------

            log.write("\n")
            log.write("=" * 78 + "\n")
            log.write(
                f"Finished: {timestamp()}\n"
            )
            log.write(
                f"Exit code: {return_code}\n"
            )
            log.write("=" * 78 + "\n")

        print()

        if return_code == 0:

            print("=" * 78)
            print(
                f"[OK] STAGE {stage['number']:02d} COMPLETED"
            )
            print("=" * 78)

            return True

        else:

            print("=" * 78)
            print(
                f"[ERROR] STAGE {stage['number']:02d} FAILED"
            )
            print(
                f"Exit code: {return_code}"
            )
            print(
                f"Log file: {log_file}"
            )
            print("=" * 78)

            return False

    except KeyboardInterrupt:

        print()
        print(
            f"[STOPPED] Stage {stage['number']:02d} "
            "interrupted by user."
        )

        try:
            process.terminate()
        except Exception:
            pass

        return False

    except Exception as e:

        print()
        print(
            f"[ERROR] Stage {stage['number']:02d} "
            f"could not be executed."
        )

        print(
            f"Reason: {e}"
        )

        return False


# =============================================================================
# PIPELINE SUMMARY
# =============================================================================

def print_pipeline_summary():
    """Print final pipeline summary."""

    print()
    print("=" * 78)
    print("PIPELINE COMPLETED")
    print("=" * 78)

    print()
    print("All stages completed successfully.")

    print()
    print("Results directory:")
    print(
        f"    {OUTPUT_DIR}"
    )

    print()
    print("Log directory:")
    print(
        f"    {LOG_DIR}"
    )

    print()
    print("Finished:", timestamp())

    print()
    print("=" * 78)


# =============================================================================
# PIPELINE FAILURE
# =============================================================================

def print_pipeline_failure(stage):
    """Print pipeline failure message."""

    print()
    print("=" * 78)
    print("PIPELINE STOPPED")
    print("=" * 78)

    print()
    print(
        f"Failed stage : {stage['number']:02d}"
    )

    print(
        f"Stage name   : {stage['name']}"
    )

    print(
        f"Script       : {stage['script']}"
    )

    print()

    if not stage.get("interactive", False):

        log_file = get_log_file(stage)

        print(
            f"Check log file:"
        )

        print(
            f"    {log_file}"
        )

    print()
    print(
        "The remaining stages were not executed."
    )

    print()
    print("=" * 78)


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():

    args = parse_args()

    # -------------------------------------------------------------------------
    # Utility modes that don't run the pipeline
    # -------------------------------------------------------------------------

    if args.reset:
        reset_state()
        sys.exit(0)

    if args.list:
        print_stage_list()
        sys.exit(0)

    # -------------------------------------------------------------------------
    # Main header
    # -------------------------------------------------------------------------

    print_header()

    # -------------------------------------------------------------------------
    # Create directories
    # -------------------------------------------------------------------------

    ensure_directories()

    # -------------------------------------------------------------------------
    # Check active runtime environment
    # -------------------------------------------------------------------------

    if not check_runtime_environment():

        sys.exit(1)

    # -------------------------------------------------------------------------
    # Check pipeline scripts
    # -------------------------------------------------------------------------

    if not check_pipeline_scripts():

        sys.exit(1)

    # -------------------------------------------------------------------------
    # Check underlying dependencies
    # -------------------------------------------------------------------------

    if not check_underlying_scripts():

        sys.exit(1)

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # Check MEME tools
    # -------------------------------------------------------------------------

    if not check_meme_tools():

        sys.exit(1)

    # -------------------------------------------------------------------------
    # Check / set up motif reference databases (JASPAR + UniProbe)
    # -------------------------------------------------------------------------

    if not check_motif_databases():

        sys.exit(1)

    # -------------------------------------------------------------------------
    # Check OrthoFinder (informational only -- see function docstring for
    # why this never blocks the pipeline the way the checks above do)
    # -------------------------------------------------------------------------

    check_orthofinder_tools()

    # -------------------------------------------------------------------------
    # Pre-flight complete
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("PRE-FLIGHT CHECKS COMPLETED")
    print("=" * 78)

    print()
    print(
        "[OK] Active Python runtime and required tools available"
    )

    print(
        "[OK] Pipeline scripts available"
    )

    print(
        "[OK] Underlying scripts available"
    )

    print(
        "[OK] MEME Suite tools available on PATH"
    )

    print(
        "[OK] MEME Suite tools available"
    )

    print(
        "[OK] Motif reference databases available"
    )

    print()

    # =========================================================================
    # DETERMINE START POINT (--resume / --from / normal full run)
    # =========================================================================

    available_numbers = sorted(s["number"] for s in STAGES)

    if args.from_stage is not None:

        start_index = resolve_start_stage(args.from_stage)

    elif args.resume:

        last = get_last_completed_stage()

        if last is None:

            print("[NOTE] No saved progress found -- starting from the beginning.")

            start_index = 0

        elif last >= available_numbers[-1]:

            print()
            print("=" * 78)
            print("Nothing left to resume -- all stages already completed.")
            print("=" * 78)

            sys.exit(0)

        else:

            print(f"[NOTE] Resuming after Stage {last:02d} (last completed).")

            start_index = resolve_start_stage(last + 1)

    else:

        start_index = 0

    start_stage = STAGES[start_index]

    # =========================================================================
    # RUN STAGE 01 (only when actually starting there)
    # =========================================================================

    if start_index == 0:

        stage_01 = STAGES[0]

        success = run_stage_01_interactive(
            stage_01
        )

        if not success:

            print_pipeline_failure(
                stage_01
            )

            sys.exit(1)

        save_state(stage_01["number"])

        print()
        print("=" * 78)
        print("STAGE 01 FINISHED")
        print("=" * 78)

        print()
        print(
            "Starting remaining pipeline stages automatically..."
        )

        print()

        remaining_stages = STAGES[1:]

    else:

        print("=" * 78)
        print("SKIPPING STAGE 01 (interactive retrieval)")
        print("=" * 78)

        print()
        print(
            "Assuming Stage 01 was already completed in a previous run. "
            "Starting directly from:"
        )
        print(
            f"    Stage {start_stage['number']:02d}: {start_stage['name']}"
        )
        print()

        remaining_stages = STAGES[start_index:]

    # =========================================================================
    # NO-PROTEIN ADJUSTMENT
    #
    # Reads outputs/.pipeline_data_manifest.json (written by Stage 01 just
    # above, or by a previous run if Stage 01 was skipped). If no real
    # protein data was found, this enforces NO_PROTEIN_MIN_STAGE as a floor
    # (unless the person explicitly passed --from, which is treated as a
    # deliberate override) and unconditionally drops any stage marked
    # "requires_protein": True in STAGES, wherever it sits.
    # =========================================================================

    manifest = load_data_manifest()

    if not manifest.get("has_real_protein", True):

        print()
        print("=" * 78)
        print("NO REAL PROTEIN DATA DETECTED")
        print("=" * 78)
        print(f"(see {DATA_MANIFEST_PATH})")
        print()

        floor_index = resolve_start_stage(NO_PROTEIN_MIN_STAGE)
        floor_number = STAGES[floor_index]["number"]

        if args.from_stage is None and floor_index > start_index:

            below_floor = [s for s in remaining_stages if s["number"] < floor_number]
            for s in below_floor:
                print(f"[SKIP] Stage {s['number']:02d}: {s['name']}  "
                      f"(below NO_PROTEIN_MIN_STAGE={NO_PROTEIN_MIN_STAGE:02d})")

            remaining_stages = [s for s in remaining_stages if s["number"] >= floor_number]

            print(f"[NOTE] Enforcing minimum start stage {NO_PROTEIN_MIN_STAGE:02d} "
                  f"because no protein data was detected.")

        elif args.from_stage is not None:

            print(f"[NOTE] --from {args.from_stage} was given explicitly -- honoring it as-is, "
                  f"not enforcing the Stage {NO_PROTEIN_MIN_STAGE:02d} floor. Still skipping any "
                  f"stage marked requires_protein below, regardless.")

        protein_only = [s for s in remaining_stages if s.get("requires_protein")]
        for s in protein_only:
            print(f"[SKIP] Stage {s['number']:02d}: {s['name']}  (requires_protein=True)")

        remaining_stages = [s for s in remaining_stages if not s.get("requires_protein")]

        if not remaining_stages:
            print()
            print("[NOTE] Nothing left to run after no-protein adjustment.")

        print()

    # =========================================================================
    # RUN REMAINING STAGES AUTOMATICALLY
    # =========================================================================

    for stage in remaining_stages:

        success = run_automatic_stage(
            stage
        )

        if not success:

            print_pipeline_failure(
                stage
            )

            print()
            print("To resume after fixing this, run either:")
            print(f"    python {Path(sys.argv[0]).name} --from {stage['number']}")
            print(f"    python {Path(sys.argv[0]).name} --resume")
            print()

            sys.exit(1)

        save_state(stage["number"])

    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================

    print_pipeline_summary()


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print()
        print()
        print("=" * 78)
        print("PIPELINE INTERRUPTED BY USER")
        print("=" * 78)

        sys.exit(130)

    except Exception as e:

        print()
        print("=" * 78)
        print("UNEXPECTED PIPELINE ERROR")
        print("=" * 78)

        print()
        print(
            f"Error: {e}"
        )

        print()
        print("=" * 78)

        sys.exit(1)
