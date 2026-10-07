#!/usr/bin/env python3
"""
Environment check for the PromEvo pipeline.

Run it with:   promevo-check
Exit code 0 = everything the pipeline needs is present, 1 = something is missing.
"""
import importlib
import os
import shutil
import sys
from pathlib import Path

PY_MODULES = [            # (import name, conda package name)
    ("numpy", "numpy"), ("pandas", "pandas"), ("scipy", "scipy"),
    ("matplotlib", "matplotlib"), ("seaborn", "seaborn"), ("Bio", "biopython"),
    ("openpyxl", "openpyxl"), ("sklearn", "scikit-learn"), ("requests", "requests"),
]

TOOLS = [                 # (label, accepted executable names)
    ("datasets (ncbi-datasets-cli)", ["datasets"]),
    ("mafft", ["mafft"]),
    ("meme", ["meme"]),
    ("tomtom", ["tomtom"]),
    ("ceqlogo", ["ceqlogo"]),
    ("orthofinder", ["orthofinder", "orthofinder.py"]),
    ("diamond", ["diamond"]),
    ("mcl", ["mcl"]),
    ("FastTree", ["FastTree", "fasttree", "FastTreeMP"]),
]

HERE = Path(__file__).resolve().parent


def main() -> int:
    problems = 0
    print("=" * 62)
    print(" PromEvo : environment check")
    print("=" * 62)
    print(f"Python {sys.version.split()[0]}  ({sys.executable})")

    print("\nPython packages")
    for mod, pkg in PY_MODULES:
        try:
            m = importlib.import_module(mod)
            print(f"  [OK]      {pkg:<14} {getattr(m, '__version__', '')}")
        except Exception as exc:
            problems += 1
            print(f"  [MISSING] {pkg:<14} ({type(exc).__name__}: {exc})")

    print("\nExternal programs on PATH")
    for label, names in TOOLS:
        found = next((shutil.which(n) for n in names if shutil.which(n)), None)
        if found:
            print(f"  [OK]      {label:<30} {found}")
        else:
            problems += 1
            print(f"  [MISSING] {label:<30} (looked for: {', '.join(names)})")

    print("\nInstalled pipeline files")
    scripts = HERE / "scripts"
    n_py = len(list(scripts.glob("*.py"))) if scripts.is_dir() else 0
    if n_py:
        print(f"  [OK]      {n_py} Python scripts in {scripts}")
    else:
        problems += 1
        print(f"  [MISSING] no .py scripts in {scripts}")
    masters = sorted(scripts.glob("master_script_*.py")) if scripts.is_dir() else []
    if masters:
        print(f"  [OK]      master script: {masters[-1].name}")
    else:
        problems += 1
        print("  [MISSING] master_script_*.py")
    # Motif databases are NOT bundled in the package: Stage 00 downloads them on the first run.
    # Same lookup order as 00_setup_motif_databases.py / 07c:
    #   1. $COMPARATIVE_PIPELINE_DATABASES   2. $XDG_DATA_HOME/comparative-gene-promoter/databases
    #   3. ~/.local/share/comparative-gene-promoter/databases
    candidates = []
    if os.environ.get("COMPARATIVE_PIPELINE_DATABASES"):
        candidates.append(Path(os.environ["COMPARATIVE_PIPELINE_DATABASES"]).expanduser())
    if os.environ.get("XDG_DATA_HOME"):
        candidates.append(Path(os.environ["XDG_DATA_HOME"]) / "comparative-gene-promoter" / "databases")
    candidates.append(Path.home() / ".local" / "share" / "comparative-gene-promoter" / "databases")
    needed = ("JASPAR2026_CORE.meme", "UniProbe_Combined.meme")
    found = next((d for d in candidates if all((d / n).is_file() and (d / n).stat().st_size > 0
                                               for n in needed)), None)
    if found:
        print(f"  [OK]      motif databases (JASPAR + UniProbe) in {found}")
    else:
        # Not a failure: the package does not ship them and the conda-build test has no network.
        print("  [INFO]    motif databases not downloaded yet -- Stage 00 fetches JASPAR + UniProbe "
              "on the first pipeline run (internet needed once)")
        print(f"            they will be stored in: {candidates[0]}")

    print("\n" + "=" * 62)
    print(" RESULT:", "all checks passed" if problems == 0 else f"{problems} problem(s) found")
    print("=" * 62)
    return 0 if problems == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
