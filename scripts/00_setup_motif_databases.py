#!/usr/bin/env python3
"""Stage 00 — one-time JASPAR + UniProbe motif database setup."""

import os
import sys
import glob
import zipfile
import subprocess
import urllib.request
from pathlib import Path

# =============================================================================
# DYNAMIC PATH RESOLUTION FOR CONDA INSTALLED PACKAGES
# =============================================================================
# 1. Determine script location
_THIS_DIR = Path(__file__).resolve().parent

# 2. Check if COMPARATIVE_PIPELINE_SCRIPTS or COMPARATIVE_PIPELINE_WORKDIR is set
ENV_SCRIPTS = os.environ.get("COMPARATIVE_PIPELINE_SCRIPTS")

if ENV_SCRIPTS:
    SCRIPTS_DIR = Path(ENV_SCRIPTS).resolve()
else:
    SCRIPTS_DIR = _THIS_DIR

# 3. Resolve parent package directory and databases directory
if SCRIPTS_DIR.name == "scripts":
    PACKAGE_DIR = SCRIPTS_DIR.parent
else:
    PACKAGE_DIR = SCRIPTS_DIR

DATABASES_DIR = PACKAGE_DIR / "databases"
DATABASES_DIR.mkdir(parents=True, exist_ok=True)

# 4. Ensure a symlink exists from scripts/databases -> databases/
# This prevents Stage 09 or other stage scripts from failing if they look in scripts/databases
LEGACY_SCRIPTS_DB_DIR = PACKAGE_DIR / "scripts" / "databases"
if PACKAGE_DIR.joinpath("scripts").exists():
    try:
        LEGACY_SCRIPTS_DB_DIR.parent.mkdir(parents=True, exist_ok=True)
        if not LEGACY_SCRIPTS_DB_DIR.exists() and not LEGACY_SCRIPTS_DB_DIR.is_symlink():
            LEGACY_SCRIPTS_DB_DIR.symlink_to(DATABASES_DIR, target_is_directory=True)
    except Exception as exc:
        # Non-critical fallback warning
        pass

# =============================================================================
# DATABASE URLS AND FILE PATHS
# =============================================================================
JASPAR_FASTA_URL = (
    "https://jaspar.elixir.no/download/data/2026/CORE/"
    "JASPAR2026_CORE_non-redundant_pfms_meme.txt"
)
JASPAR_DB_PATH = os.path.join(DATABASES_DIR, "JASPAR2026_CORE.meme")

UNIPROBE_ZIP_URL = (
    "https://thebrain.bwh.harvard.edu/uniprobe/downloads/All/All_PWMs.zip"
)
UNIPROBE_ZIP_PATH = os.path.join(DATABASES_DIR, "UniProbe_All_PWMs.zip")
UNIPROBE_EXTRACT_DIR = os.path.join(DATABASES_DIR, "uniprobe_pwms")
UNIPROBE_DB_PATH = os.path.join(DATABASES_DIR, "UniProbe_Combined.meme")


def file_is_valid(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def remove_incomplete_file(path):
    if os.path.exists(path):
        try:
            os.remove(path)
        except OSError as exc:
            print(f"[!] Could not remove incomplete file {path}: {exc}")


def download_jaspar_db():
    if file_is_valid(JASPAR_DB_PATH):
        print(f"[+] JASPAR DB already present at {JASPAR_DB_PATH} -- skipping download.")
        return True

    print("[+] Downloading JASPAR 2026 CORE motifs database...")
    try:
        subprocess.run(["wget", "-O", JASPAR_DB_PATH, JASPAR_FASTA_URL], check=True)
        if file_is_valid(JASPAR_DB_PATH):
            print(f"[+] JASPAR DB saved to {JASPAR_DB_PATH}")
            return True
        remove_incomplete_file(JASPAR_DB_PATH)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[!] Wget failed or was not found. Trying Python urllib fallback...")

    try:
        request = urllib.request.Request(
            JASPAR_FASTA_URL,
            headers={"User-Agent": "Mozilla/5.0 (compatible; ComparativeGenePromoter/0.1)"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            data = response.read()
        if not data:
            raise RuntimeError("Downloaded JASPAR file is empty.")
        with open(JASPAR_DB_PATH, "wb") as out_file:
            out_file.write(data)
        if file_is_valid(JASPAR_DB_PATH):
            print(f"[+] JASPAR DB saved to {JASPAR_DB_PATH}")
            return True
    except Exception as exc:
        print(f"[!] Failed to download JASPAR via urllib: {exc}")

    remove_incomplete_file(JASPAR_DB_PATH)
    return False


def _normalize_row(vals):
    total = sum(vals)
    if total <= 0:
        return [0.25, 0.25, 0.25, 0.25]
    return [value / total for value in vals]


def parse_pwm_file_any_orientation(file_path):
    raw_rows = []
    try:
        with open(file_path, "r", encoding="utf-8", errors="ignore") as infile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith(":"):
                    continue
                parts = line.replace(":", " ").split()
                if not parts:
                    continue
                if parts[0].upper() in ("A", "C", "G", "T"):
                    parts = parts[1:]
                nums = []
                for part in parts:
                    try:
                        nums.append(float(part))
                    except ValueError:
                        continue
                if len(nums) >= 2:
                    raw_rows.append(nums)
    except Exception:
        return None

    if not raw_rows:
        return None

    if (
        len(raw_rows) == 4
        and len(set(len(row) for row in raw_rows)) == 1
        and len(raw_rows[0]) != 4
    ):
        width = len(raw_rows[0])
        return [
            _normalize_row([raw_rows[base][position] for base in range(4)])
            for position in range(width)
        ]

    if all(len(row) == 4 for row in raw_rows):
        return [_normalize_row(row) for row in raw_rows]

    return None


def convert_uniprobe_pwms_to_meme(pwm_dir, output_meme_file):
    pwm_files = []
    for extension in ("*.pwm", "*.bml", "*.txt"):
        pwm_files.extend(
            glob.glob(os.path.join(pwm_dir, "**", extension), recursive=True)
        )
    pwm_files = sorted(set(pwm_files))
    print(f"[+] Converting {len(pwm_files)} UniProbe PWM matrices into MEME format...")

    converted_count = 0
    failed_count = 0

    try:
        with open(output_meme_file, "w", encoding="utf-8") as out:
            out.write("MEME version 4\n\n")
            out.write("ALPHABET= ACGT\n\n")
            out.write("strands: + -\n\n")

            for filepath in pwm_files:
                motif_id = os.path.splitext(os.path.basename(filepath))[0]
                try:
                    matrix = parse_pwm_file_any_orientation(filepath)
                except Exception:
                    matrix = None

                if not matrix:
                    failed_count += 1
                    continue

                out.write(f"MOTIF {motif_id}\n")
                out.write(
                    f"letter-probability matrix: alength= 4 w= {len(matrix)} nsites= 20\n"
                )
                for row in matrix:
                    out.write(
                        f"  {row[0]:.6f}\t{row[1]:.6f}\t"
                        f"{row[2]:.6f}\t{row[3]:.6f}\n"
                    )
                out.write("\n")
                converted_count += 1
    except Exception as exc:
        print(f"[!] Failed to write UniProbe MEME database: {exc}")
        remove_incomplete_file(output_meme_file)
        return False

    print(
        f"[+] Successfully generated MEME database with {converted_count} "
        f"UniProbe motifs at {output_meme_file}"
    )
    if failed_count:
        print(f"[!] {failed_count} file(s) could not be parsed as a PWM.")

    if converted_count == 0:
        remove_incomplete_file(output_meme_file)
        return False

    return file_is_valid(output_meme_file)


def download_and_setup_uniprobe():
    if file_is_valid(UNIPROBE_DB_PATH):
        print(
            f"[+] UniProbe MEME DB already exists at {UNIPROBE_DB_PATH} "
            "-- skipping download."
        )
        return True

    if not os.path.isdir(UNIPROBE_EXTRACT_DIR):
        print("[+] Downloading UniProbe PWM database from Harvard...")
        try:
            request = urllib.request.Request(
                UNIPROBE_ZIP_URL,
                headers={"User-Agent": "Mozilla/5.0 (compatible; ComparativeGenePromoter/0.1)"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                data = response.read()
            if not data:
                raise RuntimeError("Downloaded UniProbe archive is empty.")
            with open(UNIPROBE_ZIP_PATH, "wb") as out_file:
                out_file.write(data)
            print("[+] UniProbe download completed.")
        except Exception as exc:
            print(f"[!] Failed to download UniProbe database: {exc}")
            remove_incomplete_file(UNIPROBE_ZIP_PATH)
            return False

        print("[+] Extracting UniProbe PWM archive...")
        try:
            os.makedirs(UNIPROBE_EXTRACT_DIR, exist_ok=True)
            with zipfile.ZipFile(UNIPROBE_ZIP_PATH, "r") as zip_ref:
                zip_ref.extractall(UNIPROBE_EXTRACT_DIR)
            print("[+] UniProbe archive extracted.")
        except Exception as exc:
            print(f"[!] Failed to extract UniProbe archive: {exc}")
            return False
    else:
        print(f"[+] Existing UniProbe PWM extraction directory found: {UNIPROBE_EXTRACT_DIR}")

    return convert_uniprobe_pwms_to_meme(
        UNIPROBE_EXTRACT_DIR, UNIPROBE_DB_PATH
    )


def main():
    print("=" * 78)
    print("  MOTIF REFERENCE DATABASE SETUP (one-time, shared across all runs)")
    print("=" * 78)
    print(f"[*] Script location  : {_THIS_DIR}")
    print(f"[*] Package location : {PACKAGE_DIR}")
    print(f"[*] Database location: {DATABASES_DIR}")
    print()

    download_jaspar_db()
    print()
    download_and_setup_uniprobe()

    print()
    print("=" * 78)

    if file_is_valid(JASPAR_DB_PATH) and file_is_valid(UNIPROBE_DB_PATH):
        print("[✓] Motif databases ready:")
        print(f"    JASPAR   -> {JASPAR_DB_PATH}")
        print(f"    UniProbe -> {UNIPROBE_DB_PATH}")
        print("=" * 78)
        return 0

    print("[!] One or both motif databases failed to set up.")
    print(f"    JASPAR   -> {'OK' if file_is_valid(JASPAR_DB_PATH) else 'MISSING'}")
    print(f"    UniProbe -> {'OK' if file_is_valid(UNIPROBE_DB_PATH) else 'MISSING'}")
    print("=" * 78)
    return 1


if __name__ == "__main__":
    sys.exit(main())