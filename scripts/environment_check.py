#!/usr/bin/env python3
"""
Environment check for the comparative-gene-promoter pipeline.

Run it with:   comparative-gene-promoter-check
Exit code 0 = everything the pipeline needs is present, 1 = something is missing.
"""
import importlib
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
    print(" comparative-gene-promoter : environment check")
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
    # databases live next to the scripts (07c looks in <script folder>/databases)
    dbs = next((d for d in (scripts / "databases", HERE / "databases")
                if d.is_dir() and any(d.iterdir())), scripts / "databases")
    if dbs.is_dir() and any(dbs.iterdir()):
        print(f"  [OK]      motif databases in {dbs}")
    else:
        print(f"  [WARN]    no motif databases in {dbs} (07c needs them; "
              f"00_setup_motif_databases.py can create them)")

    print("\n" + "=" * 62)
    print(" RESULT:", "all checks passed" if problems == 0 else f"{problems} problem(s) found")
    print("=" * 62)
    return 0 if problems == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
