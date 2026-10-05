#!/usr/bin/env python3
"""
STAGE 03 — Build a Parallelized Neighbor-Joining (NJ) Tree with 1000 Bootstraps.

1. Calculates p-distance matrix from the MAFFT alignment.
2. Constructs a primary Neighbor-Joining (NJ) tree via Bio.Phylo.
3. Uses multiprocessing to run 1000 bootstrap replicates across all CPU threads.
4. Maps support values onto internal branches and plots a high-resolution tree.
"""

import os
import csv
import sys
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt

from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

ALIGNED_PATH = "outputs/family_aligned.fasta"
LABELS_PATH = "outputs/family_labels.csv"
OUT_TREE_IMAGE = "outputs/phylogenetic_nj_tree.png"
OUT_NEWICK_TREE = "outputs/phylogenetic_nj_tree.nwk"
NUM_BOOTSTRAPS = 10


def clean_gene_name(raw_str):
    """Extracts the full clean gene name/ID from raw headers."""
    if not raw_str:
        return ""
    parts = raw_str.split("_")
    if len(parts) >= 2 and len(parts[0]) in (6, 10) and parts[0].isalnum():
        return "_".join(parts[1:])
    return raw_str


def load_gene_labels(labels_path):
    """Loads mapping from raw FASTA header/seq_id to clean gene display name."""
    mapping = {}
    if os.path.exists(labels_path):
        with open(labels_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq_id = row.get("seq_id", "")
                gene = row.get("gene", "").strip()
                ensembl_id = row.get("ensembl_transcript_id", "").strip()

                if ensembl_id:
                    display_name = ensembl_id
                elif gene:
                    display_name = gene
                else:
                    display_name = clean_gene_name(seq_id)

                mapping[seq_id] = clean_gene_name(display_name)
    return mapping


def build_nj_tree(alignment):
    """Builds a single Neighbor-Joining tree using p-distance."""
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, method="nj")
    return constructor.build_tree(alignment)


def _single_bootstrap_worker(args):
    """
    Worker function executed in parallel across CPU cores.
    Resamples alignment columns and returns the list of leaf bipartitions.
    """
    seq_records_data, align_len, num_seqs, seed = args
    # Set unique random seed per process
    np.random.seed(seed)

    # Column resampling with replacement
    resampled_indices = np.random.choice(align_len, size=align_len, replace=True)

    resampled_records = []
    for seq_id, seq_str in seq_records_data:
        boot_seq = "".join([seq_str[i] for i in resampled_indices])
        resampled_records.append(SeqRecord(Seq(boot_seq), id=seq_id, description=""))

    resampled_align = MultipleSeqAlignment(resampled_records)
    boot_tree = build_nj_tree(resampled_align)

    # Extract bipartition sets (clades) from the replicate tree
    boot_clades = set()
    for clade in boot_tree.get_nonterminals():
        leaves = tuple(sorted([leaf.name for leaf in clade.get_terminals()]))
        if 1 < len(leaves) < num_seqs:
            boot_clades.add(leaves)

    return boot_clades


def run_parallel_bootstrap_nj(alignment, num_bootstraps=1000):
    """Performs bootstrapping concurrently using Python's multiprocessing Pool."""
    print(f"[*] Building primary Neighbor-Joining (NJ) tree...")
    main_tree = build_nj_tree(alignment)

    for clade in main_tree.get_nonterminals():
        clade.bootstrap_count = 0

    align_len = alignment.get_alignment_length()
    num_seqs = len(alignment)

    # Pre-extract sequence strings to avoid IPC serialization overhead
    seq_records_data = [(rec.id, str(rec.seq)) for rec in alignment]

    # Determine CPU worker count (Leave 1 thread free for system stability)
    num_cpus = max(1, mp.cpu_count() - 1)
    print(f"[*] Dispatching {num_bootstraps} bootstrap replicates across {num_cpus} CPU threads...")

    # Prepare argument tuples for each worker
    worker_args = [
        (seq_records_data, align_len, num_seqs, np.random.randint(0, 1_000_000_000))
        for _ in range(num_bootstraps)
    ]

    # Execute workers in parallel
    completed = 0
    with mp.Pool(processes=num_cpus) as pool:
        for boot_clades in pool.imap_unordered(_single_bootstrap_worker, worker_args, chunksize=10):
            completed += 1
            if completed % 100 == 0 or completed == num_bootstraps:
                print(f"  ➜ Progress: {completed}/{num_bootstraps} replicates done...")

            # Map bipartitions back to main tree nodes
            for clade in main_tree.get_nonterminals():
                main_leaves = tuple(sorted([leaf.name for leaf in clade.get_terminals()]))
                if main_leaves in boot_clades:
                    clade.bootstrap_count += 1

    # Calculate percentage support
    for clade in main_tree.get_nonterminals():
        confidence = (clade.bootstrap_count / num_bootstraps) * 100
        clade.confidence = round(confidence, 1)

    return main_tree


def plot_and_save_nj_tree(tree, label_mapping, out_path):
    """Renders a publication-ready NJ tree with node bootstrap values."""
    num_taxa = len(tree.get_terminals())

    # Relabel terminal leaf nodes with clean gene names
    for leaf in tree.get_terminals():
        leaf.name = label_mapping.get(leaf.name, clean_gene_name(leaf.name))

    # Dynamically scale figure height according to taxon count
    fig_height = max(12, num_taxa * 0.3)
    font_size = 5.5 if num_taxa > 150 else 7.5

    fig = plt.figure(figsize=(16, fig_height), dpi=300)
    ax = fig.add_subplot(1, 1, 1)

    # Set global leaf font size for Matplotlib cleanly
    plt.rc('font', size=font_size)

    # Plot tree via Bio.Phylo (Fixed: removed leaf_font_size kwarg)
    Phylo.draw(
        tree,
        do_show=False,
        axes=ax,
        show_confidence=True,  # Displays bootstrap support on branches
        label_func=lambda x: x.name if x.is_terminal() else ""
    )

    # Styling adjustments
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # FIX: Spine objects don't have set_title(); the axis label belongs on
    # the x-axis itself (branch-length axis), not on the spine.
    ax.set_xlabel("Substitution Distance (p-distance)")

    ax.set_title(
        f"Neighbor-Joining (NJ) Phylogenetic Tree ({num_taxa} Taxa)\n"
        f"Node values reflect Bootstrap Support (% over {NUM_BOOTSTRAPS} iterations)",
        fontsize=12,
        pad=15
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[✓] High-resolution NJ tree saved -> {out_path}")


if __name__ == "__main__":
    # Multiprocessing start method for Linux/WSL compatibility
    mp.set_start_method("spawn", force=True)

    label_mapping = load_gene_labels(LABELS_PATH)

    if os.path.exists(ALIGNED_PATH) and os.path.getsize(ALIGNED_PATH) > 0:
        print(f"[*] Loading MAFFT alignment from: {ALIGNED_PATH}")
        alignment = AlignIO.read(ALIGNED_PATH, "fasta")

        # Build parallelized NJ tree with Bootstrapping
        nj_tree = run_parallel_bootstrap_nj(alignment, num_bootstraps=NUM_BOOTSTRAPS)

        # Save tree file in Newick format
        Phylo.write(nj_tree, OUT_NEWICK_TREE, "newick")
        print(f"[✓] Newick tree structure saved -> {OUT_NEWICK_TREE}")

        # Render and save plot
        plot_and_save_nj_tree(nj_tree, label_mapping, OUT_TREE_IMAGE)

    else:
        sys.exit(f"[!] Error: Alignment file '{ALIGNED_PATH}' not found or empty. Please check Stage 02 output.")

    print("\nNext step: Run '04_conservation_and_protein_motifs.py'")