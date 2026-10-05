#!/usr/bin/env python3
"""
STAGE 00 — One-time JASPAR + UniProbe motif database setup.

This stage downloads and prepares the external motif databases required
by downstream motif-analysis stages.

IMPORTANT
---------
The motif databases are NOT packaged inside the Conda package.

Instead, they are downloaded at runtime into a user-writable directory.

Database location priority:

1. COMPARATIVE_PIPELINE_DATABASES environment variable
2. XDG_DATA_HOME/comparative-gene-promoter/databases
3. ~/.local/share/comparative-gene-promoter/databases

This makes the pipeline portable across:
    - Conda environments
    - Linux
    - WSL
    - HPC systems
    - different users
    - different installation locations
"""

import os
import sys
import glob
import zipfile
import subprocess
import urllib.request
from pathlib import Path


# =============================================================================
# DATABASE DIRECTORY RESOLUTION
# =============================================================================

def resolve_database_directory():
    """
    Determine a persistent, user-writable directory for motif databases.

    Priority:
        1. COMPARATIVE_PIPELINE_DATABASES
        2. XDG_DATA_HOME
        3. ~/.local/share
    """

    # -------------------------------------------------------------------------
    # Option 1: User explicitly specifies database location
    # -------------------------------------------------------------------------
    env_db = os.environ.get("COMPARATIVE_PIPELINE_DATABASES")

    if env_db:
        db_dir = Path(env_db).expanduser().resolve()
        db_dir.mkdir(parents=True, exist_ok=True)
        return db_dir

    # -------------------------------------------------------------------------
    # Option 2: XDG_DATA_HOME
    # -------------------------------------------------------------------------
    xdg_data_home = os.environ.get("XDG_DATA_HOME")

    if xdg_data_home:
        db_dir = (
            Path(xdg_data_home).expanduser()
            / "comparative-gene-promoter"
            / "databases"
        )
    else:
        # ---------------------------------------------------------------------
        # Option 3: Standard user data directory
        # ---------------------------------------------------------------------
        db_dir = (
            Path.home()
            / ".local"
            / "share"
            / "comparative-gene-promoter"
            / "databases"
        )

    db_dir.mkdir(parents=True, exist_ok=True)

    return db_dir


DATABASES_DIR = resolve_database_directory()


# =============================================================================
# DATABASE URLS
# =============================================================================

JASPAR_FASTA_URL = (
    "https://jaspar.elixir.no/download/data/2026/CORE/"
    "JASPAR2026_CORE_non-redundant_pfms_meme.txt"
)

UNIPROBE_ZIP_URL = (
    "https://thebrain.bwh.harvard.edu/uniprobe/downloads/All/All_PWMs.zip"
)


# =============================================================================
# DATABASE FILE PATHS
# =============================================================================

JASPAR_DB_PATH = DATABASES_DIR / "JASPAR2026_CORE.meme"

UNIPROBE_ZIP_PATH = DATABASES_DIR / "UniProbe_All_PWMs.zip"

UNIPROBE_EXTRACT_DIR = DATABASES_DIR / "uniprobe_pwms"

UNIPROBE_DB_PATH = DATABASES_DIR / "UniProbe_Combined.meme"


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def file_is_valid(path):
    """
    Return True if path exists, is a regular file, and is non-empty.
    """
    path = Path(path)

    return path.is_file() and path.stat().st_size > 0


def remove_incomplete_file(path):
    """
    Remove an incomplete or corrupted downloaded file.
    """

    path = Path(path)

    if path.exists():

        try:
            path.unlink()

        except OSError as exc:

            print(
                f"[!] Could not remove incomplete file "
                f"{path}: {exc}"
            )


# =============================================================================
# JASPAR DOWNLOAD
# =============================================================================

def download_jaspar_db():
    """
    Download the JASPAR 2026 CORE non-redundant MEME database.

    wget is attempted first because it provides robust HTTP handling.

    Python urllib is used as a fallback if wget is unavailable or fails.
    """

    if file_is_valid(JASPAR_DB_PATH):

        print(
            "[+] JASPAR database already present:"
        )

        print(
            f"    {JASPAR_DB_PATH}"
        )

        print(
            "[+] Skipping JASPAR download."
        )

        return True

    print(
        "[+] Downloading JASPAR 2026 CORE motif database..."
    )

    # -------------------------------------------------------------------------
    # Attempt wget
    # -------------------------------------------------------------------------

    try:

        subprocess.run(
            [
                "wget",
                "-O",
                str(JASPAR_DB_PATH),
                JASPAR_FASTA_URL,
            ],
            check=True,
        )

        if file_is_valid(JASPAR_DB_PATH):

            print(
                "[+] JASPAR database successfully downloaded:"
            )

            print(
                f"    {JASPAR_DB_PATH}"
            )

            return True

        remove_incomplete_file(JASPAR_DB_PATH)

    except (
        subprocess.CalledProcessError,
        FileNotFoundError,
    ):

        print(
            "[!] wget failed or is unavailable."
        )

        print(
            "[+] Trying Python urllib fallback..."
        )

    # -------------------------------------------------------------------------
    # urllib fallback
    # -------------------------------------------------------------------------

    try:

        request = urllib.request.Request(
            JASPAR_FASTA_URL,
            headers={
                "User-Agent":
                    "Mozilla/5.0 "
                    "(compatible; ComparativeGenePromoter/0.1)"
            },
        )

        with urllib.request.urlopen(
            request,
            timeout=120,
        ) as response:

            data = response.read()

        if not data:

            raise RuntimeError(
                "Downloaded JASPAR file is empty."
            )

        with open(
            JASPAR_DB_PATH,
            "wb",
        ) as out_file:

            out_file.write(data)

        if file_is_valid(JASPAR_DB_PATH):

            print(
                "[+] JASPAR database successfully downloaded:"
            )

            print(
                f"    {JASPAR_DB_PATH}"
            )

            return True

    except Exception as exc:

        print(
            f"[!] Failed to download JASPAR using urllib: {exc}"
        )

    remove_incomplete_file(JASPAR_DB_PATH)

    return False


# =============================================================================
# PWM NORMALIZATION
# =============================================================================

def _normalize_row(vals):
    """
    Normalize four nucleotide values so that their sum equals 1.
    """

    total = sum(vals)

    if total <= 0:

        return [
            0.25,
            0.25,
            0.25,
            0.25,
        ]

    return [
        value / total
        for value in vals
    ]


# =============================================================================
# UNIPROBE PWM PARSER
# =============================================================================

def parse_pwm_file_any_orientation(file_path):
    """
    Parse UniProbe PWM files regardless of whether matrices are supplied as:

        position x nucleotide

    or:

        nucleotide x position
    """

    raw_rows = []

    try:

        with open(
            file_path,
            "r",
            encoding="utf-8",
            errors="ignore",
        ) as infile:

            for line in infile:

                line = line.strip()

                if not line:
                    continue

                if line.startswith(":"):
                    continue

                parts = line.replace(
                    ":",
                    " ",
                ).split()

                if not parts:
                    continue

                # Remove nucleotide label if present
                if parts[0].upper() in (
                    "A",
                    "C",
                    "G",
                    "T",
                ):

                    parts = parts[1:]

                nums = []

                for part in parts:

                    try:

                        nums.append(
                            float(part)
                        )

                    except ValueError:

                        continue

                if len(nums) >= 2:

                    raw_rows.append(nums)

    except Exception:

        return None

    if not raw_rows:

        return None

    # -------------------------------------------------------------------------
    # Orientation 1:
    #
    # A
    # C
    # G
    # T
    #
    # where each row contains positions
    # -------------------------------------------------------------------------

    if (
        len(raw_rows) == 4
        and len(
            set(
                len(row)
                for row in raw_rows
            )
        ) == 1
        and len(raw_rows[0]) != 4
    ):

        width = len(raw_rows[0])

        return [

            _normalize_row(
                [
                    raw_rows[base][position]
                    for base in range(4)
                ]
            )

            for position in range(width)
        ]

    # -------------------------------------------------------------------------
    # Orientation 2:
    #
    # position x A,C,G,T
    # -------------------------------------------------------------------------

    if all(
        len(row) == 4
        for row in raw_rows
    ):

        return [
            _normalize_row(row)
            for row in raw_rows
        ]

    return None


# =============================================================================
# UNIPROBE → MEME CONVERSION
# =============================================================================

def convert_uniprobe_pwms_to_meme(
    pwm_dir,
    output_meme_file,
):
    """
    Convert UniProbe PWM matrices into a single MEME-format database.
    """

    pwm_dir = Path(pwm_dir)
    output_meme_file = Path(output_meme_file)

    pwm_files = []

    for extension in (
        "*.pwm",
        "*.bml",
        "*.txt",
    ):

        pwm_files.extend(
            glob.glob(
                str(pwm_dir / "**" / extension),
                recursive=True,
            )
        )

    pwm_files = sorted(
        set(pwm_files)
    )

    print(
        f"[+] Converting {len(pwm_files)} "
        "UniProbe PWM matrices into MEME format..."
    )

    converted_count = 0
    failed_count = 0

    try:

        with open(
            output_meme_file,
            "w",
            encoding="utf-8",
        ) as out:

            out.write(
                "MEME version 4\n\n"
            )

            out.write(
                "ALPHABET= ACGT\n\n"
            )

            out.write(
                "strands: + -\n\n"
            )

            for filepath in pwm_files:

                motif_id = Path(
                    filepath
                ).stem

                try:

                    matrix = (
                        parse_pwm_file_any_orientation(
                            filepath
                        )
                    )

                except Exception:

                    matrix = None

                if not matrix:

                    failed_count += 1

                    continue

                out.write(
                    f"MOTIF {motif_id}\n"
                )

                out.write(
                    "letter-probability matrix: "
                    f"alength= 4 "
                    f"w= {len(matrix)} "
                    "nsites= 20\n"
                )

                for row in matrix:

                    out.write(
                        f"  {row[0]:.6f}\t"
                        f"{row[1]:.6f}\t"
                        f"{row[2]:.6f}\t"
                        f"{row[3]:.6f}\n"
                    )

                out.write("\n")

                converted_count += 1

    except Exception as exc:

        print(
            "[!] Failed to write UniProbe MEME database:"
            f" {exc}"
        )

        remove_incomplete_file(
            output_meme_file
        )

        return False

    print(
        "[+] Successfully generated UniProbe MEME database "
        f"with {converted_count} motifs:"
    )

    print(
        f"    {output_meme_file}"
    )

    if failed_count:

        print(
            f"[!] {failed_count} file(s) "
            "could not be parsed as PWM matrices."
        )

    if converted_count == 0:

        remove_incomplete_file(
            output_meme_file
        )

        return False

    return file_is_valid(
        output_meme_file
    )


# =============================================================================
# UNIPROBE DOWNLOAD + SETUP
# =============================================================================

def download_and_setup_uniprobe():
    """
    Download UniProbe PWM archive, extract it, and convert all PWMs
    into a single MEME-format database.
    """

    if file_is_valid(
        UNIPROBE_DB_PATH
    ):

        print(
            "[+] UniProbe MEME database already exists:"
        )

        print(
            f"    {UNIPROBE_DB_PATH}"
        )

        print(
            "[+] Skipping UniProbe download."
        )

        return True

    # -------------------------------------------------------------------------
    # Download archive if necessary
    # -------------------------------------------------------------------------

    if not UNIPROBE_EXTRACT_DIR.is_dir():

        print(
            "[+] Downloading UniProbe PWM database..."
        )

        try:

            request = urllib.request.Request(
                UNIPROBE_ZIP_URL,
                headers={
                    "User-Agent":
                        "Mozilla/5.0 "
                        "(compatible; ComparativeGenePromoter/0.1)"
                },
            )

            with urllib.request.urlopen(
                request,
                timeout=120,
            ) as response:

                data = response.read()

            if not data:

                raise RuntimeError(
                    "Downloaded UniProbe archive is empty."
                )

            with open(
                UNIPROBE_ZIP_PATH,
                "wb",
            ) as out_file:

                out_file.write(data)

            print(
                "[+] UniProbe archive downloaded:"
            )

            print(
                f"    {UNIPROBE_ZIP_PATH}"
            )

        except Exception as exc:

            print(
                "[!] Failed to download UniProbe database:"
                f" {exc}"
            )

            remove_incomplete_file(
                UNIPROBE_ZIP_PATH
            )

            return False

        # ---------------------------------------------------------------------
        # Extract archive
        # ---------------------------------------------------------------------

        print(
            "[+] Extracting UniProbe PWM archive..."
        )

        try:

            UNIPROBE_EXTRACT_DIR.mkdir(
                parents=True,
                exist_ok=True,
            )

            with zipfile.ZipFile(
                UNIPROBE_ZIP_PATH,
                "r",
            ) as zip_ref:

                zip_ref.extractall(
                    UNIPROBE_EXTRACT_DIR
                )

            print(
                "[+] UniProbe archive extracted."
            )

        except Exception as exc:

            print(
                "[!] Failed to extract UniProbe archive:"
                f" {exc}"
            )

            return False

    else:

        print(
            "[+] Existing UniProbe PWM extraction directory found:"
        )

        print(
            f"    {UNIPROBE_EXTRACT_DIR}"
        )

    # -------------------------------------------------------------------------
    # Convert PWM matrices
    # -------------------------------------------------------------------------

    return convert_uniprobe_pwms_to_meme(
        UNIPROBE_EXTRACT_DIR,
        UNIPROBE_DB_PATH,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print()
    print("=" * 78)
    print(
        "  MOTIF REFERENCE DATABASE SETUP"
    )
    print(
        "  JASPAR 2026 CORE + UniProbe"
    )
    print("=" * 78)

    print(
        f"[*] Database location:"
    )

    print(
        f"    {DATABASES_DIR}"
    )

    print()

    print(
        "[*] JASPAR:"
    )

    print(
        f"    {JASPAR_DB_PATH}"
    )

    print()

    print(
        "[*] UniProbe:"
    )

    print(
        f"    {UNIPROBE_DB_PATH}"
    )

    print()
    print("-" * 78)

    # -------------------------------------------------------------------------
    # JASPAR
    # -------------------------------------------------------------------------

    jaspar_ok = download_jaspar_db()

    print()
    print("-" * 78)

    # -------------------------------------------------------------------------
    # UniProbe
    # -------------------------------------------------------------------------

    uniprobe_ok = download_and_setup_uniprobe()

    print()
    print("=" * 78)

    # -------------------------------------------------------------------------
    # Final validation
    # -------------------------------------------------------------------------

    jaspar_valid = file_is_valid(
        JASPAR_DB_PATH
    )

    uniprobe_valid = file_is_valid(
        UNIPROBE_DB_PATH
    )

    if (
        jaspar_valid
        and uniprobe_valid
    ):

        print(
            "[✓] Motif databases are ready."
        )

        print()
        print(
            f"    JASPAR   : {JASPAR_DB_PATH}"
        )

        print(
            f"    UniProbe : {UNIPROBE_DB_PATH}"
        )

        print()
        print(
            "[✓] Stage 00 completed successfully."
        )

        print("=" * 78)

        return 0

    print(
        "[!] One or both motif databases failed to set up."
    )

    print()

    print(
        "    JASPAR   : "
        + (
            "OK"
            if jaspar_valid
            else "MISSING"
        )
    )

    print(
        "    UniProbe : "
        + (
            "OK"
            if uniprobe_valid
            else "MISSING"
        )
    )

    print()
    print(
        "[!] Please check the download messages above."
    )

    print("=" * 78)

    return 1


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":

    sys.exit(
        main()
    )
