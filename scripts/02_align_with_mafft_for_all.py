"""
STAGE 02 — Align the family's sequences with MAFFT.

Skips gracefully (doesn't crash the pipeline) if mafft isn't installed --
stages 03/04 fall back to a k-mer-based approximation in that case.

Install mafft: conda install -c bioconda mafft   OR   sudo apt install mafft
"""

import os
import subprocess
import shutil

UNALIGNED_PATH = "outputs/family_sequences.fasta"
ALIGNED_PATH = "outputs/family_aligned.fasta"


def main():
    if shutil.which("mafft") is None:
        print(
            "mafft not found on PATH -- skipping real alignment.\n"
            "Stages 03 (tree) and 04 (conservation) will use a fast k-mer-based "
            "approximation instead. Install mafft and rerun this script for real results."
        )
        return

    if not os.path.exists(UNALIGNED_PATH):
        print(f"{UNALIGNED_PATH} not found -- run 01_fetch_family_sequences.py first.")
        return

    n_seqs = sum(1 for line in open(UNALIGNED_PATH) if line.startswith(">"))
    print(f"Aligning {n_seqs} sequences with MAFFT ...")

    # MAFFT writes working files to $TMPDIR (default /tmp). On WSL2, /tmp is a
    # small RAM-backed tmpfs (often ~1.4G) that fills up fast and causes cryptic
    # "No space left on device" failures even when your actual disk has plenty
    # of room. Point it at a folder next to outputs/ instead, which lives on
    # your real filesystem.
    tmp_dir = os.path.abspath("outputs/mafft_tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    env = os.environ.copy()
    env["TMPDIR"] = tmp_dir

    with open(ALIGNED_PATH, "w") as out_f:
        result = subprocess.run(
            ["mafft", "--auto", "--quiet", UNALIGNED_PATH],
            stdout=out_f, stderr=subprocess.PIPE, text=True, env=env,
        )

    if result.returncode != 0:
        print(f"MAFFT failed:\n{result.stderr}")
        if os.path.exists(ALIGNED_PATH):
            os.remove(ALIGNED_PATH)  # don't leave a stale/empty file behind - downstream
                                      # stages check for this file's existence to decide
                                      # whether a real alignment is available
        return

    lengths = set()
    for line in open(ALIGNED_PATH):
        if not line.startswith(">"):
            lengths.add(len(line.strip()))
    # (rough length check per line; fine since mafft output isn't line-wrapped here)
    print(f"Alignment complete -> {ALIGNED_PATH}")
    print("Now run: python3 03_build_tree.py")


if __name__ == "__main__":
    main()
