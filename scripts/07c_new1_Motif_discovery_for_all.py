#!/usr/bin/env python3
"""
STAGE 07c — MEME Motif Discovery + TOMTOM Database Search (per-organism)

This stage:

1. Reads outputs/promoter_sequences.fasta.
2. Uses outputs/family_labels.csv to determine the organism for each
   promoter sequence.
3. Runs MEME de novo motif discovery separately for each organism.
4. Runs TOMTOM against:
       - JASPAR2026 CORE
       - UniProbe
5. Performs an additional genus-level MEME/TOMTOM analysis.
6. Produces combined TOMTOM summary files.
7. Produces a complete MEME run-status report.

Motif databases are NOT packaged with the Conda installation.

They are expected at:

    ~/.local/share/comparative-gene-promoter/databases/

or, if specified:

    COMPARATIVE_PIPELINE_DATABASES=/path/to/databases

Stage 00 is responsible for downloading and preparing the databases.
Stage 07c only validates and reads them.
"""


import os
import re
import sys
import glob
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# =============================================================================
# PROGRAM INFORMATION
# =============================================================================

PROGRAM_NAME = "Comparative Gene Promoter"
PROGRAM_VERSION = "0.1.0"


# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FASTA = "outputs/promoter_sequences.fasta"

FAMILY_LABELS_CSV = "outputs/family_labels.csv"

BASE_OUT_DIR = "outputs"


# =============================================================================
# SPECIES-LEVEL OUTPUTS
# =============================================================================

MEME_OUT_DIR = os.path.join(
    BASE_OUT_DIR,
    "meme_output"
)

TOMTOM_OUT_DIR = os.path.join(
    BASE_OUT_DIR,
    "tomtom_output"
)

PER_ORG_FASTA_DIR = os.path.join(
    BASE_OUT_DIR,
    "promoters_by_organism_for_meme"
)


# =============================================================================
# GENUS-LEVEL OUTPUTS
# =============================================================================

MEME_OUT_DIR_GENUS = os.path.join(
    BASE_OUT_DIR,
    "meme_output_genus"
)

TOMTOM_OUT_DIR_GENUS = os.path.join(
    BASE_OUT_DIR,
    "tomtom_output_genus"
)

PER_GENUS_FASTA_DIR = os.path.join(
    BASE_OUT_DIR,
    "promoters_by_genus_for_meme"
)

GENUS_SUMMARY_CSV = os.path.join(
    BASE_OUT_DIR,
    "combined_tomtom_summary_genus.csv"
)


# =============================================================================
# MEME SETTINGS
# =============================================================================

MEME_NUM_MOTIFS = 5

MEME_MINW = 6

MEME_MAXW = 15

MEME_USE_REVCOMP = True

MEME_MAXSIZE_FLOOR = 10_000_000

MIN_SEQS_FOR_MEME = 2

EXPECTED_PROMOTER_LEN = 1000


# =============================================================================
# STATUS REPORT
# =============================================================================

RUN_STATUS_CSV = os.path.join(
    BASE_OUT_DIR,
    "meme_run_status.csv"
)


# =============================================================================
# DYNAMIC SCRIPT PATH RESOLUTION
# =============================================================================

_THIS_DIR = Path(__file__).resolve().parent


ENV_SCRIPTS = os.environ.get(
    "COMPARATIVE_PIPELINE_SCRIPTS"
)


if ENV_SCRIPTS:

    SCRIPTS_DIR = (
        Path(ENV_SCRIPTS)
        .expanduser()
        .resolve()
    )

else:

    SCRIPTS_DIR = _THIS_DIR


# =============================================================================
# PACKAGE DIRECTORY
# =============================================================================

if SCRIPTS_DIR.name == "scripts":

    PACKAGE_DIR = SCRIPTS_DIR.parent

else:

    PACKAGE_DIR = SCRIPTS_DIR


# =============================================================================
# DATABASE DIRECTORY RESOLUTION
# =============================================================================
#
# The database location must be consistent with Stage 00.
#
# Priority:
#
# 1. COMPARATIVE_PIPELINE_DATABASES
# 2. XDG_DATA_HOME/comparative-gene-promoter/databases
# 3. ~/.local/share/comparative-gene-promoter/databases
#
# For backwards compatibility, existing package-level database directories
# are also checked AFTER the user cache location.
#
# IMPORTANT:
# 07c does NOT create the database directory.
# Stage 00 creates/downloads the databases.
#


def resolve_database_directory():

    # -------------------------------------------------------------------------
    # 1. Explicit user-defined database location
    # -------------------------------------------------------------------------

    env_database_dir = os.environ.get(
        "COMPARATIVE_PIPELINE_DATABASES"
    )

    if env_database_dir:

        return (
            Path(env_database_dir)
            .expanduser()
            .resolve()
        )


    # -------------------------------------------------------------------------
    # 2. XDG cache location
    # -------------------------------------------------------------------------

    xdg_cache_home = os.environ.get(
        "XDG_DATA_HOME"
    )

    if xdg_cache_home:

        xdg_database_dir = (
            Path(xdg_cache_home)
            .expanduser()
            .resolve()
            / "comparative-gene-promoter"
            / "databases"
        )

    else:

        xdg_database_dir = (
            Path.home()
            / ".cache"
            / "comparative-gene-promoter"
            / "databases"
        )


    # -------------------------------------------------------------------------
    # 3. Preferred user cache directory
    # -------------------------------------------------------------------------

    if xdg_database_dir.exists():

        return xdg_database_dir.resolve()


    # -------------------------------------------------------------------------
    # 4. Legacy package-level database directory
    #
    # This is retained for compatibility with older installations.
    # -------------------------------------------------------------------------

    legacy_candidates = [

        PACKAGE_DIR / "databases",

        SCRIPTS_DIR / "databases",

        Path(
            os.environ.get(
                "CONDA_PREFIX",
                ""
            )
        )
        / "share"
        / "comparative-gene-promoter"
        / "databases",
    ]


    for candidate in legacy_candidates:

        if candidate.exists() and candidate.is_dir():

            return candidate.resolve()


    # -------------------------------------------------------------------------
    # 5. Return preferred cache location even if it does not yet exist.
    #
    # Stage 07c will report that the databases are missing.
    # It will NOT create the directory.
    # -------------------------------------------------------------------------

    return xdg_database_dir.resolve()


DATABASES_DIR = resolve_database_directory()


# =============================================================================
# DATABASE FILE PATHS
# =============================================================================

JASPAR_DB_PATH = str(
    DATABASES_DIR /
    "JASPAR2026_CORE.meme"
)

UNIPROBE_DB_PATH = str(
    DATABASES_DIR /
    "UniProbe_Combined.meme"
)


# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

os.makedirs(
    BASE_OUT_DIR,
    exist_ok=True
)

os.makedirs(
    MEME_OUT_DIR,
    exist_ok=True
)

os.makedirs(
    TOMTOM_OUT_DIR,
    exist_ok=True
)

os.makedirs(
    PER_ORG_FASTA_DIR,
    exist_ok=True
)

os.makedirs(
    MEME_OUT_DIR_GENUS,
    exist_ok=True
)

os.makedirs(
    TOMTOM_OUT_DIR_GENUS,
    exist_ok=True
)

os.makedirs(
    PER_GENUS_FASTA_DIR,
    exist_ok=True
)


# =============================================================================
# DATABASE VALIDATION
# =============================================================================

def check_reference_databases():
    """
    Check that JASPAR and UniProbe databases exist and are non-empty.

    Stage 07c does not download databases.
    Database installation is handled by Stage 00.
    """

    missing = []

    databases = [
        (
            "JASPAR2026_CORE",
            JASPAR_DB_PATH
        ),
        (
            "UniProbe_Combined",
            UNIPROBE_DB_PATH
        ),
    ]


    for label, path in databases:

        if (
            not os.path.isfile(path)
            or os.path.getsize(path) == 0
        ):

            missing.append(
                (label, path)
            )


    # -------------------------------------------------------------------------
    # Databases available
    # -------------------------------------------------------------------------

    if not missing:

        print(
            "[+] Reference databases found at:"
        )

        print(
            f"    {DATABASES_DIR}"
        )

        print(
            f"      JASPAR   -> {JASPAR_DB_PATH}"
        )

        print(
            f"      UniProbe -> {UNIPROBE_DB_PATH}"
        )

        return True


    # -------------------------------------------------------------------------
    # Databases missing
    # -------------------------------------------------------------------------

    print()
    print(
        "[!] Missing reference database(s)."
    )

    print(
        f"    Expected database directory:"
    )

    print(
        f"    {DATABASES_DIR}"
    )

    for label, path in missing:

        print(
            f"      {label} -> NOT FOUND ({path})"
        )


    print()
    print(
        "[!] The motif databases must be created by Stage 00."
    )

    print()
    print(
        "[+] Run Stage 00 before running Stage 07c:"
    )

    print(
        f"    python "
        f"{os.path.join(_THIS_DIR, '00_setup_motif_databases.py')}"
    )

    print()
    print(
        "[+] Alternatively, specify a custom database location:"
    )

    print(
        "    export COMPARATIVE_PIPELINE_DATABASES=/path/to/databases"
    )

    print()

    return False


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def sanitize_filename(name):

    safe = re.sub(
        r"[^A-Za-z0-9]+",
        "_",
        (name or "").strip()
    )

    return (
        safe.strip("_")
        or "Unknown_organism"
    )


def extract_genus(organism):

    organism = (
        organism or ""
    ).strip()

    if (
        not organism
        or organism.lower() == "nan"
    ):

        return "Unknown_genus"

    return organism.split()[0]


# =============================================================================
# ORGANISM MAPPING
# =============================================================================

def load_organism_map(labels_path):

    if not os.path.exists(labels_path):

        print(
            f"[!] Family labels file not found: {labels_path}"
        )

        return {}


    try:

        df = pd.read_csv(
            labels_path
        )

    except Exception as exc:

        print(
            f"[!] Could not read {labels_path}: {exc}"
        )

        return {}


    lower_cols = {
        str(c).lower(): c
        for c in df.columns
    }


    seq_id_col = lower_cols.get(
        "seq_id"
    )


    organism_col = next(
        (
            lower_cols[c]
            for c in (
                "organism",
                "species",
                "organism_name"
            )
            if c in lower_cols
        ),
        None
    )


    if (
        seq_id_col is None
        or organism_col is None
    ):

        print(
            f"[!] Could not identify seq_id and organism "
            f"columns in {labels_path}"
        )

        return {}


    return dict(
        zip(
            df[seq_id_col].astype(str),
            df[organism_col].astype(str)
        )
    )


# =============================================================================
# FASTA READING & GROUPING
# =============================================================================

def read_and_validate_fasta(fasta_file):

    if not os.path.exists(fasta_file):

        print(
            f"[!] Input FASTA file '{fasta_file}' not found."
        )

        return []


    try:

        records = list(
            SeqIO.parse(
                fasta_file,
                "fasta"
            )
        )

    except Exception as exc:

        print(
            f"[!] Could not parse FASTA: {exc}"
        )

        return []


    print(
        f"[+] Successfully loaded "
        f"{len(records)} sequences from {fasta_file}."
    )

    return records


def group_records_by_organism(
    records,
    organism_map
):

    groups = {}


    for rec in records:

        if organism_map:

            organism = organism_map.get(
                rec.id,
                "Unknown_organism"
            )

        else:

            organism = "All_organisms"


        groups.setdefault(
            organism,
            []
        ).append(rec)


    return groups


def group_records_by_genus(
    records,
    organism_map
):

    groups = {}


    for rec in records:

        if organism_map:

            organism = organism_map.get(
                rec.id,
                ""
            )

            genus = extract_genus(
                organism
            )

        else:

            genus = "All_organisms"


        groups.setdefault(
            genus,
            []
        ).append(rec)


    return groups


def write_organism_fasta(
    records,
    out_path
):

    with open(
        out_path,
        "w"
    ) as outfile:

        for rec in records:

            outfile.write(
                f">{rec.id}\n"
            )

            outfile.write(
                f"{str(rec.seq)}\n"
            )


def count_meme_motifs(
    meme_file_path
):

    if not os.path.exists(
        meme_file_path
    ):

        return 0


    count = 0


    try:

        with open(
            meme_file_path,
            "r",
            encoding="utf-8",
            errors="ignore"
        ) as infile:

            for line in infile:

                if line.startswith(
                    "MOTIF"
                ):

                    count += 1


    except Exception:

        return 0


    return count


# =============================================================================
# RUN MEME
# =============================================================================

def run_meme(
    fasta_file,
    out_dir,
    num_motifs=MEME_NUM_MOTIFS,
    minw=MEME_MINW,
    maxw=MEME_MAXW
):

    try:

        total_chars = sum(
            len(rec.seq)
            for rec in SeqIO.parse(
                fasta_file,
                "fasta"
            )
        )

    except Exception as exc:

        return (
            None,
            f"meme_failed: could not read FASTA: {exc}"
        )


    maxsize = max(
        MEME_MAXSIZE_FLOOR,
        int(total_chars * 1.1)
    )


    meme_cmd = [
        "meme",
        fasta_file,
        "-dna",
        "-mod",
        "zoops",
        "-nmotifs",
        str(num_motifs),
        "-minw",
        str(minw),
        "-maxw",
        str(maxw),
        "-maxsize",
        str(maxsize),
    ]


    if MEME_USE_REVCOMP:

        meme_cmd.append(
            "-revcomp"
        )


    meme_cmd += [
        "-oc",
        out_dir
    ]


    # -------------------------------------------------------------------------
    # Remove stale MEME output files
    # -------------------------------------------------------------------------

    for stale in (
        "meme.txt",
        "meme.xml",
        "meme.html"
    ):

        stale_path = os.path.join(
            out_dir,
            stale
        )

        if os.path.exists(
            stale_path
        ):

            try:

                os.remove(
                    stale_path
                )

            except OSError:

                pass


    os.makedirs(
        out_dir,
        exist_ok=True
    )


    print(
        f"[+] Executing MEME on {fasta_file}"
    )


    try:

        res = subprocess.run(
            meme_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

    except FileNotFoundError:

        print(
            "[!] MEME executable not found in PATH."
        )

        return (
            None,
            "meme_failed: executable not found"
        )


    meme_txt = os.path.join(
        out_dir,
        "meme.txt"
    )


    if (
        res.returncode != 0
        or not os.path.exists(meme_txt)
    ):

        tail = (
            res.stderr or ""
        ).strip()[-1000:]


        print(
            f"[!] MEME failed "
            f"(exit code {res.returncode})."
        )


        return (
            None,
            f"meme_failed: exit "
            f"{res.returncode}: {tail}"
        )


    n_found = count_meme_motifs(
        meme_txt
    )


    print(
        f"[+] MEME execution complete "
        f"({n_found} motif(s) found)."
    )


    if n_found == 0:

        return (
            None,
            "meme_zero_motifs"
        )


    return (
        meme_txt,
        "ok"
    )


# =============================================================================
# RUN TOMTOM
# =============================================================================

def run_tomtom_search(
    meme_file,
    db_path,
    db_name,
    out_dir
):

    if (
        not os.path.exists(db_path)
        or os.path.getsize(db_path) == 0
    ):

        print(
            f"[!] Database missing or empty: {db_path}"
        )

        return None


    if count_meme_motifs(
        db_path
    ) == 0:

        print(
            f"[!] {db_name} contains 0 motifs."
        )

        return None


    os.makedirs(
        out_dir,
        exist_ok=True
    )


    stale_tsv = os.path.join(
        out_dir,
        "tomtom.tsv"
    )


    if os.path.exists(
        stale_tsv
    ):

        try:

            os.remove(
                stale_tsv
            )

        except OSError:

            pass


    tomtom_cmd = [
        "tomtom",
        "-oc",
        out_dir,
        "-thresh",
        "10.0",
        "-evalue",
        "-min-overlap",
        "4",
        meme_file,
        db_path,
    ]


    print(
        f"[+] Running TOMTOM against {db_name}..."
    )


    try:

        res = subprocess.run(
            tomtom_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

    except FileNotFoundError:

        print(
            "[!] TOMTOM executable not found in PATH."
        )

        return None


    if res.returncode != 0:

        tail = (
            res.stderr or ""
        ).strip()[-1000:]


        print(
            f"[!] TOMTOM failed against "
            f"{db_name} "
            f"(exit code {res.returncode})."
        )


        if tail:

            print(
                f"[!] TOMTOM error: {tail}"
            )


    tsv_path = os.path.join(
        out_dir,
        "tomtom.tsv"
    )


    return (
        tsv_path
        if os.path.exists(tsv_path)
        else None
    )


# =============================================================================
# TOMTOM COLUMN FINDER
# =============================================================================

def _find_column(
    df,
    *candidates
):

    lower_map = {
        str(c).lower(): c
        for c in df.columns
    }


    for cand in candidates:

        if cand.lower() in lower_map:

            return lower_map[
                cand.lower()
            ]


    return None


# =============================================================================
# PARSE TOMTOM TSV
# =============================================================================

def parse_tomtom_tsv(
    tomtom_tsv,
    db_source
):

    if (
        not tomtom_tsv
        or not os.path.exists(tomtom_tsv)
    ):

        return pd.DataFrame()


    try:

        df = pd.read_csv(
            tomtom_tsv,
            sep="\t",
            comment="#"
        )

    except Exception as exc:

        print(
            f"[!] Error reading TOMTOM output TSV: {exc}"
        )

        return pd.DataFrame()


    if df.empty:

        return df


    query_col = _find_column(
        df,
        "query_id",
        "Query_ID"
    )


    target_col = _find_column(
        df,
        "target_id",
        "Target_ID"
    )


    if (
        query_col is None
        or target_col is None
    ):

        return pd.DataFrame()


    df = df.dropna(
        subset=[
            query_col,
            target_col
        ]
    )


    df["Database"] = db_source


    df = df.rename(
        columns={
            query_col: "Query_ID",
            target_col: "Target_ID"
        }
    )


    return df


# =============================================================================
# MEME + TOMTOM PIPELINE
# =============================================================================

def run_meme_tomtom_pipeline(
    label,
    records,
    fasta_dir,
    meme_out_root,
    tomtom_out_root,
    tag_column
):

    slug = sanitize_filename(
        label
    )


    n_seqs = len(
        records
    )


    total_bp = sum(
        len(r.seq)
        for r in records
    )


    print(
        f"\n{'=' * 70}"
    )

    print(
        f"{tag_column.capitalize()}: {label}"
    )

    print(
        f"Sequences: {n_seqs:,}"
    )

    print(
        f"Total promoter bp: {total_bp:,}"
    )

    print(
        f"{'=' * 70}"
    )


    status = {

        "level":
        tag_column,

        "group":
        label,

        "n_sequences":
        n_seqs,

        "total_bp":
        total_bp,

        "status":
        "",

        "detail":
        "",
    }


    if n_seqs < MIN_SEQS_FOR_MEME:

        msg = (
            f"only {n_seqs} sequence(s); "
            f"MEME needs at least "
            f"{MIN_SEQS_FOR_MEME}."
        )


        print(
            f"[!] {tag_column.capitalize()} "
            f"'{label}': {msg}"
        )


        status.update(
            status=
            "skipped_too_few_sequences",

            detail=
            msg
        )


        return (
            [],
            status
        )


    fasta_path = os.path.join(
        fasta_dir,
        f"{slug}.fasta"
    )


    write_organism_fasta(
        records,
        fasta_path
    )


    meme_dir = os.path.join(
        meme_out_root,
        slug
    )


    meme_file, detail = run_meme(
        fasta_path,
        meme_dir
    )


    if not meme_file:

        status.update(
            status=
            detail.split(":")[0],

            detail=
            detail
        )


        return (
            [],
            status
        )


    tomtom_dir = os.path.join(
        tomtom_out_root,
        slug
    )


    combined_results = []


    # =========================================================================
    # JASPAR
    # =========================================================================

    jaspar_tsv = run_tomtom_search(

        meme_file,

        JASPAR_DB_PATH,

        "JASPAR2026_CORE",

        os.path.join(
            tomtom_dir,
            "JASPAR2026_CORE"
        ),
    )


    if jaspar_tsv:

        df_j = parse_tomtom_tsv(
            jaspar_tsv,
            "JASPAR2026_CORE"
        )


        if not df_j.empty:

            df_j[tag_column] = label

            combined_results.append(
                df_j
            )


    # =========================================================================
    # UNIPROBE
    # =========================================================================

    uniprobe_tsv = run_tomtom_search(

        meme_file,

        UNIPROBE_DB_PATH,

        "UniProbe",

        os.path.join(
            tomtom_dir,
            "UniProbe"
        ),
    )


    if uniprobe_tsv:

        df_u = parse_tomtom_tsv(
            uniprobe_tsv,
            "UniProbe"
        )


        if not df_u.empty:

            df_u[tag_column] = label

            combined_results.append(
                df_u
            )


    total_hits = sum(
        len(d)
        for d in combined_results
    )


    status.update(

        status="ok",

        detail=
        f"{total_hits} TOMTOM hit(s)"
    )


    return (
        combined_results,
        status
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print(
        "\n" + "=" * 78
    )

    print(
        "STAGE 07c — MEME MOTIF DISCOVERY + TOMTOM"
    )

    print(
        "=" * 78
    )

    print(
        f"[*] Program version   : {PROGRAM_VERSION}"
    )

    print(
        f"[*] Script directory  : {_THIS_DIR}"
    )

    print(
        f"[*] Database directory: {DATABASES_DIR}"
    )

    print(
        f"[*] JASPAR database   : {JASPAR_DB_PATH}"
    )

    print(
        f"[*] UniProbe database : {UNIPROBE_DB_PATH}"
    )

    print(
        f"[*] Working directory : {os.getcwd()}"
    )

    print(
        "=" * 78
    )

    print()


    # =========================================================================
    # INPUT FASTA
    # =========================================================================

    records = read_and_validate_fasta(
        INPUT_FASTA
    )


    if not records:

        return 1


    # =========================================================================
    # DATABASE CHECK
    # =========================================================================

    if not check_reference_databases():

        print()

        print(
            "[ERROR] Stage 07c cannot continue because "
            "motif databases are missing."
        )

        return 1


    # =========================================================================
    # ORGANISM MAP
    # =========================================================================

    organism_map = load_organism_map(
        FAMILY_LABELS_CSV
    )


    # =========================================================================
    # GROUP BY SPECIES / ORGANISM
    # =========================================================================

    species_groups = (
        group_records_by_organism(
            records,
            organism_map
        )
    )


    # =========================================================================
    # GROUP BY GENUS
    # =========================================================================

    genus_groups = (
        group_records_by_genus(
            records,
            organism_map
        )
    )


    # =========================================================================
    # STATUS STORAGE
    # =========================================================================

    all_statuses = []


    # =========================================================================
    # SPECIES-LEVEL ANALYSIS
    # =========================================================================

    species_results = []


    print()
    print(
        "=" * 78
    )

    print(
        "SPECIES / ORGANISM-LEVEL MEME + TOMTOM ANALYSIS"
    )

    print(
        "=" * 78
    )


    for label, recs in species_groups.items():

        res, stat = (
            run_meme_tomtom_pipeline(

                label,

                recs,

                PER_ORG_FASTA_DIR,

                MEME_OUT_DIR,

                TOMTOM_OUT_DIR,

                "organism",
            )
        )


        species_results.extend(
            res
        )


        all_statuses.append(
            stat
        )


    # =========================================================================
    # COMBINED SPECIES TOMTOM SUMMARY
    # =========================================================================

    if species_results:

        combined_sp = pd.concat(
            species_results,
            ignore_index=True
        )


        combined_sp.to_csv(

            os.path.join(
                BASE_OUT_DIR,
                "combined_tomtom_summary.csv"
            ),

            index=False
        )


        print()
        print(
            "[+] Species-level TOMTOM summary written:"
        )

        print(
            f"    {os.path.join(BASE_OUT_DIR, 'combined_tomtom_summary.csv')}"
        )


    # =========================================================================
    # GENUS-LEVEL ANALYSIS
    # =========================================================================

    genus_results = []


    print()
    print(
        "=" * 78
    )

    print(
        "GENUS-LEVEL MEME + TOMTOM ANALYSIS"
    )

    print(
        "=" * 78
    )


    for label, recs in genus_groups.items():

        res, stat = (
            run_meme_tomtom_pipeline(

                label,

                recs,

                PER_GENUS_FASTA_DIR,

                MEME_OUT_DIR_GENUS,

                TOMTOM_OUT_DIR_GENUS,

                "genus",
            )
        )


        genus_results.extend(
            res
        )


        all_statuses.append(
            stat
        )


    # =========================================================================
    # COMBINED GENUS TOMTOM SUMMARY
    # =========================================================================

    if genus_results:

        combined_gn = pd.concat(
            genus_results,
            ignore_index=True
        )


        combined_gn.to_csv(
            GENUS_SUMMARY_CSV,
            index=False
        )


        print()
        print(
            "[+] Genus-level TOMTOM summary written:"
        )

        print(
            f"    {GENUS_SUMMARY_CSV}"
        )


    # =========================================================================
    # RUN STATUS
    # =========================================================================

    pd.DataFrame(
        all_statuses
    ).to_csv(
        RUN_STATUS_CSV,
        index=False
    )


    # =========================================================================
    # FINAL MESSAGE
    # =========================================================================

    print()
    print(
        "=" * 78
    )

    print(
        "STAGE 07c COMPLETE"
    )

    print(
        "=" * 78
    )

    print(
        f"[+] Status report:"
    )

    print(
        f"    {RUN_STATUS_CSV}"
    )

    print(
        f"[+] Database directory:"
    )

    print(
        f"    {DATABASES_DIR}"
    )

    print(
        "=" * 78
    )


    return 0


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":

    sys.exit(
        main()
    )