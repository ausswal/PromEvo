#!/usr/bin/env python3
"""
STAGE 01 — Combined Gene Family (protein) & Promoter Fetcher, Dual-Mode,
Interactive, Multi-Organism.

Fetches PROTEIN sequences and PROMOTER sequences for the same gene family
in one run, keeping them in two clearly separate output namespaces that are
joined only through shared seq_id / organism keys -- exactly what a
downstream comparative analysis (e.g. Step 23's species x feature matrix)
needs: one consistent identity per gene, independent data per sequence type.

  PROTEIN outputs (unchanged from 01c/01d):
    outputs/family_sequences.fasta          (master protein FASTA)
    outputs/by_organism/<Organism>.fasta    (per-organism protein FASTAs)
    outputs/family_labels.csv               (seq_id, family, length, organism,
                                              gene, domain_trimmed,
                                              ensembl_transcript_id, ncbi_gene_id)
    outputs/domain_architecture.csv

  PROMOTER outputs (always written by every mode):
    outputs/promoter_sequences.fasta
    outputs/promoter_fetch_log.csv          (seq_id, gene_id, status, organism,
                                              match_method -- same columns
                                              regardless of which mode/method
                                              produced them)
    outputs/promoters_by_organism/<Organism>/promoter_sequences.fasta
    outputs/promoters_by_organism/<Organism>/promoter_fetch_log.csv

How promoters are obtained depends on the mode:
  Mode 1 (local genome + GFF):
    promoters are extracted directly from the same local GFF pass that
    finds the genes, then split into outputs/promoters_by_organism/ for
    consistency with the other mode -- see split_local_promoters_by_organism().
    Every GFF file found under the directory you point it at (searched
    recursively, so both "one subfolder per organism" and "every genome's
    files sitting flat in one folder" work identically) is paired with its
    genome FASTA and treated as its own organism/strain -- see
    discover_gff_fasta_pairs(). The organism name is read directly from the
    genome FASTA's own header when it carries one (the VEuPathDB/ToxoDB
    convention: ">seqid | organism=Genus_species_strain | ..."), not
    guessed from a folder name, so a single flat directory containing many
    strains' files (e.g. many "ToxoDB-68_Tgondii<Strain>.gff" /
    "..._Genome.fasta" pairs) is correctly split into that many distinct
    organisms rather than lumped together.
  Mode 2 (NCBI protein online):
    proteins are fetched first; promoters are then obtained by downloading
    each distinct organism's reference genome + GFF3 via NCBI Datasets and
    matching genes locally (GeneID -> gene symbol -> family-keyword, in that
    order of preference) -- see fetch_and_write_promoters_via_genome().
    Requires the `datasets` CLI on PATH.

CORRECTNESS FIX vs. the earlier separate 01c/01d + 05 scripts: promoter
FASTA headers now always use the exact same seq_id as family_labels.csv
(built once, by make_seq_id(), and reused everywhere an ID is constructed)
-- previously, Mode 1's promoter FASTA used a bare gene_id while
family_labels.csv used "{gene_id}_{gene_id}", so promoters silently
couldn't be joined back to their organism for local-genome-derived data.
"""

import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from collections import defaultdict
from Bio import Entrez, SeqIO

from pipeline_config import PROMOTER_LEN   # the ONE shared promoter length (1000 bp)

# NCBI requires a contact email for Entrez, and rate-limits to ~3 req/sec
NCBI_EMAIL = "your_email@example.com"
NCBI_REQUEST_DELAY = 0.34
Entrez.email = NCBI_EMAIL

# Promoter geometry for the genome-download path (Mode 2). Mode 1's local
# GFF path uses its own promoter_upstream prompt instead (upstream-only,
# matching that mode's existing UX).
PROMOTER_UPSTREAM_BP = PROMOTER_LEN
PROMOTER_DOWNSTREAM_BP = 0
GENOME_ROOT = "outputs/genomes_by_organism"


def make_seq_id(acc, gene):
    """The one place seq_id is constructed, used identically for protein
    records, family_labels.csv, and promoter FASTA headers, so every output
    namespace joins on the exact same key."""
    return f"{acc}_{gene}".replace(" ", "_").replace(";", "").replace("/", "_")


def sanitize_organism_name(organism):
    """Alphanumeric-only organism -> folder-name sanitizer, used consistently
    for BOTH outputs/by_organism/ and outputs/promoters_by_organism/ (and
    matching Step 23's own sanitize() in 23_comparative_matrix.py) so a
    given organism always maps to the same folder name everywhere."""
    if not organism or str(organism).strip() == "" or str(organism).lower() == "nan":
        return "Unknown_organism"
    return re.sub(r"[^A-Za-z0-9]+", "_", str(organism).strip()).strip("_")



# ==============================================================================
# HELPER: INTERACTIVE MULTI-ORGANISM PROMPTER
# ==============================================================================

def get_organism_list_from_user(default_orgs=None):
    """
    Interactively prompts user for the number of organisms and collects their names sequentially.
    """
    if default_orgs is None:
        default_orgs = ["Plasmodium"]

    print("\n--- ORGANISM SELECTION ---")
    num_input = input("How many target species/organisms do you want to fetch? (Enter a number, or 'all' for no organism restriction) [default: 1]: ").strip()

    if num_input.lower() == "all":
        print("  ➜ Querying across ALL species (no organism filter applied).")
        return []

    if not num_input:
        num_species = len(default_orgs)
    elif num_input.isdigit() and int(num_input) > 0:
        num_species = int(num_input)
    else:
        print("  [!] Invalid input. Defaulting to 1 species.")
        num_species = 1

    org_list = []
    print(f"Please enter the names for {num_species} species/organisms:")
    for i in range(1, num_species + 1):
        default_val = default_orgs[i - 1] if i <= len(default_orgs) else ""
        prompt = f"  Species {i}" + (f" [default: {default_val}]" if default_val else "") + ": "
        val = input(prompt).strip()
        
        if not val and default_val:
            val = default_val
        if val:
            org_list.append(val)

    return org_list


def get_local_root_dirs_from_user(default_dir="Genomes"):
    """
    Interactively prompts for how many separate local directories to scan
    for genome FASTA + GFF files, then collects each path in turn -- same
    "how many, then loop" UX as get_organism_list_from_user() (used for
    Mode 2), just for directory paths instead of organism names.

    This is for when your genus folders (e.g. one for Plasmodium, one for
    Toxoplasma) do NOT share a common parent directory you could point at
    once. If they DO share a parent, you only need one directory here --
    discover_gff_fasta_pairs() already searches recursively (subfolders
    included), so a single parent folder containing multiple genus
    subfolders is auto-detected in one pass regardless.
    """
    print("\n--- LOCAL DATA DIRECTORIES ---")
    num_input = input(
        "How many separate local directories do you want to scan for genome "
        "FASTA + GFF files? (e.g. 2 if your Plasmodium and Toxoplasma data "
        "live in two unrelated folders; 1 if they share a common parent "
        "folder you can point at once -- every subfolder under it is "
        "searched automatically) [default: 1]: "
    ).strip()

    if not num_input:
        num_dirs = 1
    elif num_input.isdigit() and int(num_input) > 0:
        num_dirs = int(num_input)
    else:
        print("  [!] Invalid input. Defaulting to 1 directory.")
        num_dirs = 1

    root_dirs = []
    for i in range(1, num_dirs + 1):
        default_val = default_dir if (i == 1 and num_dirs == 1) else ""
        prompt = (f"  Directory {i}"
                  + (f" [default: {default_val}]" if default_val else "")
                  + ": ")
        while True:
            val = input(prompt).strip()
            if not val and default_val:
                val = default_val
            if not val:
                print("    [!] Path cannot be empty.")
                continue
            if not os.path.isdir(val):
                print(f"    [!] '{val}' is not a directory or doesn't exist -- try again.")
                continue
            root_dirs.append(val)
            break

    return root_dirs


def discover_gff_fasta_pairs(root_dir):
    """Recursively finds every .gff/.gff3 file under root_dir -- root_dir
    itself, or any subfolder -- so this works identically whether files are
    laid out one-subfolder-per-organism, or (as with a typical ToxoDB/
    VEuPathDB bulk download) every strain's genome FASTA + GFF sitting
    flat together in one directory.

    Each GFF is paired with its genome FASTA (see find_matching_fasta) and
    assigned an organism/strain name (see extract_organism_from_fasta_header,
    with guess_organism_from_filename as a fallback) -- so a single flat
    directory containing many strains' files is correctly split into that
    many distinct organisms, not lumped together under one folder name.

    Returns a list of {"gff": path, "fasta": path, "organism": name,
    "protein_fasta": path_or_None} dicts. A GFF with no matching genome
    FASTA is skipped with a printed warning, never silently dropped. A
    missing protein_fasta is NOT a warning -- most local datasets won't
    have one, and process_local_genome_and_gff() falls back to nucleotide
    gene sequences in that case."""
    gff_paths = sorted(
        glob.glob(os.path.join(root_dir, "**", "*.gff"), recursive=True) +
        glob.glob(os.path.join(root_dir, "**", "*.gff3"), recursive=True)
    )

    pairs = []
    for gff_path in gff_paths:
        fasta_path = find_matching_fasta(gff_path)
        if not fasta_path:
            print(f"  [!] No matching genome FASTA found for '{gff_path}' -- skipping.")
            continue

        organism = extract_organism_from_fasta_header(fasta_path)
        if not organism:
            organism = guess_organism_from_filename(gff_path, root_dir)

        protein_fasta = find_matching_protein_fasta(gff_path)
        if protein_fasta:
            print(f"  [i] Found companion protein FASTA for '{os.path.basename(gff_path)}': "
                  f"{os.path.basename(protein_fasta)}")

        pairs.append({"gff": gff_path, "fasta": fasta_path, "organism": organism,
                      "protein_fasta": protein_fasta})

    return pairs


def find_matching_protein_fasta(gff_path):
    """
    Finds a companion protein FASTA for one GFF file. Tries the
    VEuPathDB/ToxoDB naming convention first (e.g. "ToxoDB-68_TgondiiME49.gff"
    pairs with "ToxoDB-68_TgondiiME49_AnnotatedProteins.fasta"), then NCBI's
    ".faa" extension (never matched by the VEuPathDB-oriented checks below,
    since NCBI protein FASTAs are conventionally ".faa", not ".fasta"/".fa"),
    then falls back to any protein-extension FASTA in the same directory
    whose name contains "protein" (case-insensitive) and starts with the
    GFF's own stem, and finally -- only when unambiguous -- to the single
    protein-named FASTA sitting in the same directory regardless of its
    stem. That last fallback exists because NCBI Datasets CLI downloads
    commonly place a genome, its GFF, and its protein FASTA in one
    per-assembly folder with names that share NO common stem at all (e.g.
    "genomic.gff" + "genomic.fna" + "protein.faa").

    Returns None if nothing matches -- that's the normal, expected outcome
    for most local GFF/genome-FASTA pairs (most local datasets won't have
    a protein FASTA at all), not an error. process_local_genome_and_gff()
    falls back to extracting nucleotide gene sequences in that case, exactly
    as it always has.
    """
    protein_exts = (".fasta", ".fa", ".pep", ".faa")

    base_prefix = os.path.splitext(gff_path)[0]
    exact_candidates = [
        base_prefix + "_AnnotatedProteins.fasta",
        base_prefix + "_AnnotatedProteins.fa",
        base_prefix + "_Proteins.fasta",
        base_prefix + "_protein.fasta",
        base_prefix + ".protein.fasta",
        base_prefix + ".pep.fasta",
        base_prefix + ".pep",
        base_prefix + ".faa",
        base_prefix + "_protein.faa",
    ]
    for cand in exact_candidates:
        if os.path.exists(cand):
            return cand

    gff_stem = os.path.basename(base_prefix)
    directory = os.path.dirname(gff_path) or "."
    for fname in sorted(os.listdir(directory)):
        stem, ext = os.path.splitext(fname)
        if (ext.lower() in protein_exts
                and stem.startswith(gff_stem)
                and re.search(r"protein", stem, re.IGNORECASE)):
            return os.path.join(directory, fname)

    same_dir_protein_files = [
        fname for fname in sorted(os.listdir(directory))
        if os.path.splitext(fname)[1].lower() in protein_exts
        and re.search(r"protein", fname, re.IGNORECASE)
    ]
    if len(same_dir_protein_files) == 1:
        return os.path.join(directory, same_dir_protein_files[0])

    return None


def find_matching_fasta(gff_path):
    """Finds the genome FASTA that goes with one GFF file. Tries the exact
    ToxoDB/VEuPathDB convention first (e.g. "Foo.gff" + "Foo_Genome.fasta")
    and a few other common suffixes, then falls back to ANY FASTA file in
    the same directory whose name (minus extension) starts with the GFF's
    own name -- covers naming conventions not explicitly listed below
    without requiring an exact match."""
    base_prefix = os.path.splitext(gff_path)[0]
    exact_candidates = [
        base_prefix + "_Genome.fasta",
        base_prefix + ".fasta",
        base_prefix + ".fa",
        base_prefix + ".fna",
    ]
    for cand in exact_candidates:
        if os.path.exists(cand):
            return cand

    gff_stem = os.path.basename(base_prefix)
    directory = os.path.dirname(gff_path) or "."
    for fname in sorted(os.listdir(directory)):
        stem, ext = os.path.splitext(fname)
        if ext.lower() in (".fasta", ".fa", ".fna") and stem.startswith(gff_stem):
            return os.path.join(directory, fname)

    return None


def extract_organism_from_fasta_header(fasta_path):
    """Reads just the first header line of a genome FASTA and pulls out an
    organism name -- only the first line is read, so this is cheap even on
    a multi-gigabyte genome FASTA. Two conventions are tried, in order:

    1. VEuPathDB/ToxoDB's 'organism=' field, e.g.:
         >KQ973609 | organism=Toxoplasma_gondii_ARI | version=2016-03-14 ...

    2. NCBI RefSeq/GenBank's plain free-text description (no 'organism='
       field at all), e.g.:
         >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
       Here the organism name is everything after the accession, with the
       trailing assembly-boilerplate clause NCBI conventionally appends
       (", complete genome" / ", complete sequence" / ", chromosome N" /
       ", whole genome shotgun sequence" / etc.) trimmed off.

    Returns None (so the caller falls back to a filename-based guess) if
    neither convention yields anything."""
    try:
        with open(fasta_path, "r") as f:
            first_line = f.readline()
    except Exception:
        return None
    if not first_line.startswith(">"):
        return None

    m = re.search(r"organism=([^|]+)", first_line)
    if m:
        return m.group(1).strip().replace("_", " ")

    header = first_line[1:].strip()
    parts = header.split(None, 1)
    if len(parts) < 2:
        return None
    description = parts[1].strip()
    if not description:
        return None
    description = re.split(
        r",?\s*(?:complete genome|complete sequence|complete cds|"
        r"chromosome\b.*|plasmid\b.*|whole genome shotgun sequence.*|"
        r"partial genome.*|genome assembly.*|contig\b.*|scaffold\b.*|unplaced\b.*)",
        description, maxsplit=1, flags=re.IGNORECASE,
    )[0].strip().rstrip(",")
    return description or None


def guess_organism_from_filename(gff_path, root_dir):
    """Fallback organism name when the FASTA header has no 'organism='
    field: prefer the immediate parent subfolder name if the file isn't
    sitting directly in root_dir (the older one-subfolder-per-organism
    convention), otherwise fall back to the GFF's own filename, stripped of
    a leading numeric database-version prefix like 'ToxoDB-68_' and its
    extension."""
    parent = os.path.basename(os.path.dirname(gff_path))
    root_name = os.path.basename(os.path.normpath(root_dir))
    if parent and parent != root_name:
        return parent.replace("_", " ")

    stem = os.path.splitext(os.path.basename(gff_path))[0]
    stem = re.sub(r"^[A-Za-z]+DB-\d+_", "", stem)  # strip "ToxoDB-68_" style prefixes
    return stem.replace("_", " ")


def get_pairs_from_user(pairs):
    """Interactively choose how many of the discovered genome/GFF pairs to
    process -- same 'a number, or all' prompt style used elsewhere in this
    script, just operating on discovered file pairs instead of organism
    names or folders."""
    print("\n--- LOCAL GENOME/GFF FILES FOUND ---")
    if not pairs:
        return []

    print(f"Found {len(pairs)} genome FASTA + GFF pair(s):")
    for i, p in enumerate(pairs, 1):
        protein_marker = "  [protein available]" if p.get("protein_fasta") else ""
        print(f"    {i}. {p['organism']}  ({os.path.basename(p['gff'])}){protein_marker}")

    num_input = input("\nHow many of these do you want to process? "
                       "(Enter a number, or 'all') [default: all]: ").strip()

    if not num_input or num_input.lower() == "all":
        print(f"  ➜ Processing ALL {len(pairs)} pair(s).")
        return pairs

    if not num_input.isdigit() or int(num_input) <= 0:
        print("  [!] Invalid input. Defaulting to 'all'.")
        return pairs

    num = min(int(num_input), len(pairs))
    sel_input = input(f"  Enter {num} pair number(s) (space/comma separated), "
                       f"or press Enter to take the first {num} listed: ").strip()
    if not sel_input:
        return pairs[:num]

    chosen = []
    for tok in re.split(r"[,\s]+", sel_input):
        if tok.isdigit() and 1 <= int(tok) <= len(pairs):
            chosen.append(pairs[int(tok) - 1])
    if not chosen:
        print(f"  [!] No valid selection made. Defaulting to the first {num}.")
        return pairs[:num]
    return chosen


# ==============================================================================
# MODE 2: NCBI PROTEIN DATABASE ONLINE FETCHING
# ==============================================================================

def fetch_all_ncbi_ids_via_history(full_query, batch_size=500):
    """Runs one esearch to get the total count + history, then pages through
    ALL of it via retstart/retmax -- so nothing is capped by a fixed
    'max sequences' number. Returns the full id_list."""
    handle = Entrez.esearch(db="protein", term=full_query, retmax=0, usehistory="y")
    search_record = Entrez.read(handle)
    handle.close()
    time.sleep(NCBI_REQUEST_DELAY)

    total_count = int(search_record.get("Count", 0))
    if total_count == 0:
        return [], 0

    webenv, query_key = search_record["WebEnv"], search_record["QueryKey"]
    if total_count > 2000:
        print(f"    [i] {total_count} total matches — this is a large result set, "
              f"fetching all of it may take a while and is subject to NCBI rate limits.")

    ids = []
    for start in range(0, total_count, batch_size):
        try:
            handle = Entrez.esearch(db="protein", term=full_query, retstart=start,
                                     retmax=batch_size, webenv=webenv, query_key=query_key)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(NCBI_REQUEST_DELAY)
            ids.extend(record.get("IdList", []))
        except Exception as e:
            print(f"    [!] Failed to fetch ID batch at offset {start}: {e}")
    return ids, total_count


def fetch_ncbi_records_by_id(id_list):
    """Fetches FASTA + GenBank details and the linked GeneID for each protein
    ID.

    CORRECTNESS FIX (found 2026-09): also falls back to the CDS/gene
    feature's `/locus_tag` qualifier when `/gene` is absent, before ever
    falling back to the bare accession. Many non-reference-strain genome
    submissions (WGS/TSA-derived -- typical for the less-studied strains of
    apicomplexan parasites) were never assigned a curated short gene
    SYMBOL and were never indexed into NCBI's formal Gene database at all
    -- so both `/gene` and the elink(protein->gene) GeneID lookup below
    come up empty for them, no matter which NCBI endpoint you fetch from.
    Their `/locus_tag` (e.g. "TGARI_268850" / "TGGT1_320000"), however, is
    assigned by the submission pipeline for essentially every gene
    regardless of curation status -- and critically, that exact string is
    what the genome's own downloaded GFF3 uses as that gene's Name/ID
    attribute too. Using locus_tag as the "gene" column value here means
    match_genes_in_gff's gene_symbol fallback can actually succeed for
    these genes downstream, instead of trying to match a bare protein
    accession (e.g. "KYF50247") against the GFF, which never works."""
    records = []
    n_linked = 0
    n_locus_tag_fallback = 0

    for i, pid in enumerate(id_list, 1):
        try:
            fasta_handle = Entrez.efetch(db="protein", id=pid, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(fasta_handle, "fasta")
            fasta_handle.close()
            time.sleep(NCBI_REQUEST_DELAY)

            gb_handle = Entrez.efetch(db="protein", id=pid, rettype="gb", retmode="text")
            gb_record = SeqIO.read(gb_handle, "gb")
            gb_handle.close()
            time.sleep(NCBI_REQUEST_DELAY)
        except Exception as e:
            print(f"  [{i}/{len(id_list)}] {pid} — fetch failed: {e}")
            continue

        organism = gb_record.annotations.get("organism", "")
        relevant_feats = [f for f in gb_record.features if f.type in ("gene", "CDS")]

        gene_name = ""
        for feat in relevant_feats:
            if "gene" in feat.qualifiers:
                gene_name = feat.qualifiers["gene"][0]
                break

        used_locus_tag = False
        if not gene_name:
            for feat in relevant_feats:
                if "locus_tag" in feat.qualifiers:
                    gene_name = feat.qualifiers["locus_tag"][0]
                    used_locus_tag = True
                    break

        acc = seq_record.id.split(".")[0]
        if not gene_name:
            gene_name = acc
        else:
            n_locus_tag_fallback += int(used_locus_tag)

        gene_id = ""
        try:
            link_handle = Entrez.elink(dbfrom="protein", db="gene", id=pid)
            link_record = Entrez.read(link_handle)
            link_handle.close()
            time.sleep(NCBI_REQUEST_DELAY)
            linksets = link_record[0].get("LinkSetDb", [])
            if linksets and linksets[0].get("Link"):
                gene_id = linksets[0]["Link"][0]["Id"]
        except Exception:
            pass

        seq = str(seq_record.seq)
        records.append({
            "Entry": acc,
            "Sequence": seq,
            "Gene Names": gene_name,
            "Organism": organism,
            "_full_length": len(seq),
            "_all_domains": [("Full_Protein", 1, len(seq))],
            "_domain_trimmed": False,
            "_ensembl_transcript": "",
            "_ncbi_gene_id": gene_id
        })

        n_linked += int(bool(gene_id))
        if gene_id:
            status = f"GeneID {gene_id}"
        elif used_locus_tag:
            status = f"no GeneID -- using locus_tag '{gene_name}' for symbol matching"
        else:
            status = "no linked GeneID, no gene/locus_tag qualifier either"
        print(f"  [{i}/{len(id_list)}] {acc} ({organism}) — {status}")

    if n_locus_tag_fallback:
        print(f"\n  [i] {n_locus_tag_fallback}/{len(id_list)} record(s) had no /gene symbol and no "
              f"linked GeneID, but DID have a /locus_tag -- these are only matchable downstream "
              f"via gene_symbol against the GFF's Name/locus_tag attribute, not exact_geneid.")

    return records, n_linked


def fetch_ncbi_protein_online(family_query, org_list):
    """
    Fetches ALL matching NCBI Protein records for family_query, querying
    each organism in org_list SEPARATELY and exhaustively (via
    usehistory/retstart pagination) rather than one combined OR-query --
    a combined query lets whichever organism has more matching records
    dominate the result and silently starve the others. If org_list is
    empty ('all' was chosen), runs a single unrestricted query.
    """
    field_attempts = [
        ("Protein Name", "({q}[Protein Name])"),
        ("Title", "({q}[Title])"),
        ("All Fields", "({q})"),
    ]

    organisms_to_query = org_list if org_list else [None]
    all_records = []
    total_linked = 0

    for org in organisms_to_query:
        org_label = org if org else "ALL ORGANISMS (no restriction)"
        org_clause = f'"{org}"[Organism]' if org else ""

        print(f"\n{'=' * 60}\nQuerying NCBI Protein for organism: {org_label}\n{'=' * 60}")

        id_list, full_query = [], ""
        for field_label, clause_template in field_attempts:
            clause = clause_template.format(q=family_query)
            query_parts = [clause]
            if org_clause:
                query_parts.append(org_clause)
            full_query = " AND ".join(query_parts)

            print(f"[➜] Trying NCBI Protein Query ({field_label}): {full_query}")
            try:
                id_list, total_count = fetch_all_ncbi_ids_via_history(full_query)
            except Exception as e:
                print(f"  [!] Query failed ({field_label}): {e}")
                id_list, total_count = [], 0

            if id_list:
                print(f"[✓] Matched via [{field_label}] — {total_count} total record(s) available.")
                if field_label == "All Fields":
                    print("  [!] Note: this is the least precise field — matches can include "
                          "off-target hits where the keyword appears for unrelated reasons. "
                          "Review the organism/name printed for each fetched record below.")
                break
            else:
                print(f"  [!] No hits via [{field_label}], broadening search...")

        if not id_list:
            print(f"  [!] No results found for organism '{org_label}' — skipping.")
            continue

        print(f"\nFetching details for {len(id_list)} record(s) from {org_label}...\n")
        org_records, n_linked = fetch_ncbi_records_by_id(id_list)
        all_records.extend(org_records)
        total_linked += n_linked

    if not all_records:
        sys.exit("✗ No results returned for any organism — check your family keyword or organism names.")

    print(f"\n[✓] Fetched {len(all_records)} total protein record(s) from NCBI across "
          f"{len(organisms_to_query)} organism quer{'y' if len(organisms_to_query) == 1 else 'ies'} "
          f"({total_linked}/{len(all_records)} with a linked GeneID for promoter matching).")
    return all_records, family_query


# ==============================================================================
# MODE 1: LOCAL GENOME & GFF (DEDUPLICATED PROMOTERS)
# ==============================================================================

def parse_gff_attributes(attr_str):
    attrs = {}
    for item in attr_str.strip().split(";"):
        if not item.strip():
            continue
        if "=" in item:
            k, v = item.strip().split("=", 1)
            attrs[k.strip()] = v.strip()
        elif " " in item:
            parts = item.strip().split(None, 1)
            if len(parts) == 2:
                attrs[parts[0].strip()] = parts[1].strip().replace('"', '')
    return attrs


def clean_gene_id(raw_id, strip_version=True):
    """Normalize an ID to base gene level (removes transcript suffixes like -t30_1).

    strip_version controls the generic trailing ".<digits>" strip (e.g.
    "PF3D7_0100100.1" -> "PF3D7_0100100", or protein accession
    "CAB95322.1" -> "CAB95322"). That strip is right for transcript/protein
    IDs but WRONG for gene-feature IDs in some organisms: Trypanosoma
    locus tags like "Tb927.1.100" would collapse to "Tb927.1", merging
    thousands of distinct genes into ~a dozen IDs. Gene-feature IDs (and
    anything resolved to one) are therefore cleaned with
    strip_version=False; protein/transcript IDs keep the default."""
    if strip_version:
        base_id = re.sub(r'(-t\d+.*|\.t\d+.*|-mRNA.*|\.\d+)$', '', raw_id, flags=re.IGNORECASE)
    else:
        base_id = re.sub(r'(-t\d+.*|\.t\d+.*|-mRNA.*)$', '', raw_id, flags=re.IGNORECASE)
    # NCBI RefSeq/GenBank GFF3 prefixes every ID by feature type --
    # "gene-b0001", "rna-NM_...", "cds-NP_..." -- so a gene's own ID and a
    # CDS's Parent= reference to it both carry this prefix. VEuPathDB IDs
    # never use this convention, so stripping it is a safe no-op there.
    base_id = re.sub(r'^(gene|rna|cds|exon)-', '', base_id, flags=re.IGNORECASE)
    return base_id.strip()


def load_protein_lookup(protein_fasta_path):
    """
    Loads a companion protein FASTA into {base_gene_id: (sequence, header_id)},
    keyed by the SAME clean_gene_id() normalization the GFF gene-ID parsing
    already uses -- VEuPathDB protein FASTA headers are usually transcript-
    level (e.g. "TGME49_268560-t26_1"), so stripping the transcript suffix
    is what lets these match up with the GFF's gene-level IDs at all.

    When multiple transcripts collapse to the same base gene ID (isoforms),
    the LONGEST protein is kept as that gene's representative sequence --
    a standard, explicit convention, not an arbitrary pick.
    """
    lookup = {}
    if not protein_fasta_path or not os.path.exists(protein_fasta_path):
        return lookup

    n_isoform_collisions = 0
    for rec in SeqIO.parse(protein_fasta_path, "fasta"):
        header_id = rec.id
        base_id = clean_gene_id(header_id)
        seq = str(rec.seq).strip()
        if not seq:
            continue
        if base_id in lookup:
            n_isoform_collisions += 1
            if len(seq) <= len(lookup[base_id][0]):
                continue  # keep the longer isoform already stored
        lookup[base_id] = (seq, header_id)

    if n_isoform_collisions:
        print(f"    [i] {n_isoform_collisions} isoform(s) collapsed to their gene's "
              f"longest protein sequence.")

    return lookup


def build_cds_annotations(gff_path):
    """Pre-scans the GFF once to recover, per GENE, two things NCBI-style
    GFF3 'gene' features don't carry themselves:

    1. protein_id -- the accession used as the header in a companion .faa
       (e.g. "CAB95322.1"). Not derivable from the gene ID by any string
       transform; the only link is the CDS line's protein_id= attribute.
    2. product -- the human-readable description (e.g. "hypothetical
       protein, unlikely"), which NCBI puts on the CDS/mRNA line, not the
       gene line, so a keyword search over gene lines alone misses it.

    The CDS's Parent= is either the gene directly (simple prokaryote-style
    gene -> CDS) or an mRNA/transcript whose own Parent= is the gene
    (eukaryote-style gene -> mRNA -> CDS, e.g. Trypanosoma). Both are
    resolved by following the Parent chain up to a gene-type feature.

    Returns {gene_id (cleaned, strip_version=False): {"protein_id": ...,
    "product": ...}}. A multi-exon CDS / multi-isoform gene keeps the first
    protein_id/product seen."""
    annotations = {}
    if not os.path.exists(gff_path):
        return annotations

    parent_of = {}      # raw feature ID -> raw parent ID (first parent)
    is_gene = set()     # raw IDs of gene-type features
    cds_rows = []       # (parent_raw, protein_id, product)
    rna_product = {}    # raw mRNA/transcript ID -> product (fallback)

    with open(gff_path, "r") as gff_handle:
        for line in gff_handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2].lower()
            attrs = parse_gff_attributes(parts[8])
            fid = attrs.get("ID")
            parent = attrs.get("Parent", "").split(",")[0].strip()
            if ftype in ("gene", "protein_coding_gene") and fid:
                is_gene.add(fid)
            elif ftype == "cds":
                if parent:
                    cds_rows.append((parent, attrs.get("protein_id"), attrs.get("product", "")))
            elif fid and parent:
                parent_of[fid] = parent
                if ftype in ("mrna", "transcript") and attrs.get("product"):
                    rna_product[fid] = attrs["product"]

    for parent, protein_id, product in cds_rows:
        node, rna_id, hops = parent, None, 0
        while node not in is_gene and node in parent_of and hops < 5:
            rna_id = rna_id or node
            node = parent_of[node]
            hops += 1
        if node not in is_gene:
            continue  # couldn't resolve to a gene
        gene_id = clean_gene_id(node, strip_version=False)
        if not product and rna_id:
            product = rna_product.get(rna_id, "")
        if gene_id and gene_id not in annotations:
            annotations[gene_id] = {"protein_id": protein_id, "product": product}
    return annotations


def load_protein_lookup_by_exact_id(protein_fasta_path):
    """Loads a companion protein FASTA into {exact_header_id: sequence},
    with NO gene/transcript-suffix normalization -- used to match NCBI
    protein accessions (e.g. "NP_414542.1") exactly as they appear in both
    the .faa header and the GFF CDS line's protein_id= attribute. This is
    deliberately separate from load_protein_lookup()'s clean_gene_id()-based
    matching, which exists for VEuPathDB's different (suffix-based) ID
    convention and would be actively wrong here -- e.g. it strips a
    trailing ".1" version number, which is part of the actual accession."""
    lookup = {}
    if not protein_fasta_path or not os.path.exists(protein_fasta_path):
        return lookup
    for rec in SeqIO.parse(protein_fasta_path, "fasta"):
        seq = str(rec.seq).strip()
        if seq:
            lookup[rec.id] = seq
    return lookup



def process_local_genome_and_gff(gff_path, fasta_path, organism_name, gene_keyword,
                                  promoter_upstream=1000, protein_fasta_path=None):
    """Scans ONE genome's GFF for gene features matching gene_keyword, and
    extracts both the gene sequence and its upstream promoter directly from
    the paired genome FASTA. organism_name is passed in (resolved once per
    pair by discover_gff_fasta_pairs -- normally read straight from the
    genome FASTA's own 'organism=' header field) rather than guessed here.

    If protein_fasta_path is given (see find_matching_protein_fasta), each
    gene's REAL protein sequence is used instead of its raw nucleotide gene
    sequence whenever a match is found by ID -- so local runs with a
    companion protein FASTA available produce genuinely protein-typed
    output (detect_sequence_type()/has_real_protein will correctly reflect
    this), same as Mode 2. Genes with no protein match (not every gene
    model gets one -- pseudogenes, ncRNA, incomplete annotations, etc.)
    fall back to their nucleotide sequence exactly as before; that's
    normal and expected, not an error."""
    if not os.path.exists(fasta_path):
        print(f"[!] Genome FASTA '{fasta_path}' not found -- skipping '{organism_name}'.")
        return [], []

    print(f"\n[*] Scanning '{os.path.basename(gff_path)}' ({organism_name}) for keyword: '{gene_keyword}'...")

    protein_lookup = load_protein_lookup(protein_fasta_path)
    if protein_lookup:
        print(f"    [i] Loaded {len(protein_lookup)} protein sequence(s) from "
              f"'{os.path.basename(protein_fasta_path)}' for ID matching.")

    # NCBI-style companion data: protein_lookup_by_exact_id supports exact
    # accession matching (e.g. "NP_414542.1"), and cds_annotations supplies
    # both that accession and the CDS 'product' description per gene --
    # neither is derivable from the gene line alone in NCBI GFF3. Both are
    # harmless no-ops for VEuPathDB-style data (protein_lookup_by_exact_id
    # stays empty when it's not a real companion match; cds_annotations
    # simply won't have entries if there are no CDS lines with the
    # relevant attributes).
    protein_lookup_by_exact_id = load_protein_lookup_by_exact_id(protein_fasta_path)
    cds_annotations = build_cds_annotations(gff_path)

    records = []
    promoter_records = []
    seen_loci = set()
    n_protein_matched = 0

    genome_seqs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    with open(gff_path, "r") as gff_handle:
        for line in gff_handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start_str, end_str, score, strand, phase, attributes_raw = parts

            if feature_type.lower() in ["gene", "protein_coding_gene"]:
                attrs = parse_gff_attributes(attributes_raw)
                raw_gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name") or "Unknown"
                gene_id = clean_gene_id(raw_gene_id, strip_version=False)
                description = attrs.get("description") or attrs.get("Note") or attrs.get("Name") or ""

                cds_info = cds_annotations.get(gene_id, {})
                product = cds_info.get("product", "")

                full_text = f"{gene_id} {description} {product} {attributes_raw}".lower()
                if gene_keyword.lower() in full_text:
                    start, end = int(start_str), int(end_str)

                    locus_key = (seqid, start, end, strand)
                    if locus_key in seen_loci:
                        continue
                    seen_loci.add(locus_key)

                    if seqid not in genome_seqs:
                        continue

                    chr_seq = genome_seqs[seqid].seq
                    chr_len = len(chr_seq)

                    if strand == "+":
                        prom_end_coord = start - 1
                        prom_start_coord = max(0, prom_end_coord - promoter_upstream)
                        prom_seq_obj = chr_seq[prom_start_coord:prom_end_coord]
                    else:
                        prom_start_coord = end
                        prom_end_coord = min(chr_len, prom_start_coord + promoter_upstream)
                        prom_seq_forward = chr_seq[prom_start_coord:prom_end_coord]
                        prom_seq_obj = prom_seq_forward.reverse_complement()

                    if len(prom_seq_obj) > 0:
                        # Correctness fix: key promoter records by the SAME
                        # seq_id write_outputs() will construct for this row
                        # (make_seq_id(acc, gene), both = gene_id here), not
                        # by the bare gene_id -- otherwise
                        # outputs/promoter_sequences.fasta headers can't be
                        # joined back to family_labels.csv's seq_id column.
                        # NOTE: promoters are always nucleotide, regardless
                        # of whether a protein match is used for "Sequence"
                        # below -- a promoter is a regulatory DNA region,
                        # not something a protein FASTA could ever supply.
                        promoter_records.append((make_seq_id(gene_id, gene_id), str(prom_seq_obj)))

                    gene_genomic_seq = chr_seq[start-1:end]
                    if strand == "-":
                        gene_genomic_seq = gene_genomic_seq.reverse_complement()

                    protein_match = None
                    ncbi_protein_id = cds_info.get("protein_id")
                    if ncbi_protein_id and ncbi_protein_id in protein_lookup_by_exact_id:
                        protein_match = (protein_lookup_by_exact_id[ncbi_protein_id], ncbi_protein_id)
                    elif gene_id in protein_lookup:
                        protein_match = protein_lookup[gene_id]

                    if protein_match:
                        seq_out = protein_match[0]
                        is_real_protein = True
                        n_protein_matched += 1
                    else:
                        seq_out = str(gene_genomic_seq)
                        is_real_protein = False

                    records.append({
                        "Entry": gene_id,
                        "Sequence": seq_out,
                        "Gene Names": gene_id,
                        "Organism": organism_name,
                        "_full_length": len(seq_out),
                        "_domain_trimmed": False,
                        "_ensembl_transcript": gene_id,
                        "_all_domains": [("Target_Gene", 1, len(seq_out))],
                        "_ncbi_gene_id": "",
                        "_is_real_protein": is_real_protein,
                    })

    if not records:
        n_gene_lines = n_with_product = 0
        with open(gff_path, "r") as _h:
            for _l in _h:
                _p = _l.split("\t")
                if len(_p) >= 9 and _p[2].lower() in ("gene", "protein_coding_gene"):
                    n_gene_lines += 1
        n_with_product = sum(1 for v in cds_annotations.values() if v.get("product"))
        print(f"    [?] 0 matches for '{gene_keyword}'. Diagnostic: {n_gene_lines} gene line(s) in GFF, "
              f"{len(cds_annotations)} gene(s) linked to a CDS, {n_with_product} of those with a product= "
              f"description. If gene lines is 0, this GFF uses a different gene feature type; if the "
              f"linked count is 0, its CDS Parent= chain doesn't reach a gene.")

    if protein_lookup:
        print(f"    [i] {n_protein_matched}/{len(records)} gene(s) matched to a real protein "
              f"sequence by ID; the rest fall back to their nucleotide gene sequence.")

    return records, promoter_records


def process_local_genomes_multi_org(pairs, gene_keyword, promoter_upstream=1000):
    """
    Runs process_local_genome_and_gff() once per discovered
    {"gff", "fasta", "organism"} pair and aggregates the results -- lets
    Mode 1 process any number of genomes in one run (5, 10, or 'all'),
    regardless of whether they're laid out one-subfolder-per-organism or
    all sitting flat in a single directory (see discover_gff_fasta_pairs).
    """
    all_records, all_promoter_records = [], []
    n_with_protein_source = 0

    for pair in pairs:
        print(f"\n{'=' * 60}")
        print(f"Processing: {pair['organism']}  ({pair['gff']})")
        print(f"{'=' * 60}")
        records, promoter_records = process_local_genome_and_gff(
            pair["gff"], pair["fasta"], pair["organism"], gene_keyword, promoter_upstream,
            protein_fasta_path=pair.get("protein_fasta"))
        print(f"  ➜ {len(records)} matching gene(s) found in this organism.")
        n_with_protein_source += sum(1 for r in records if r.get("_is_real_protein"))
        all_records.extend(records)
        all_promoter_records.extend(promoter_records)

    print(f"\n[✓] Multi-organism extraction complete: {len(all_records)} total gene(s) "
          f"across {len(pairs)} genome(s).")
    if n_with_protein_source:
        print(f"    {n_with_protein_source}/{len(all_records)} gene(s) use a real protein "
              f"sequence (matched from a companion protein FASTA); the rest use their "
              f"nucleotide gene sequence.")
    return all_records, all_promoter_records, gene_keyword


# ==============================================================================
# PROMOTER FETCHING -- MODE 2 (genome download + local GFF matching)
# ==============================================================================
# One genome+GFF3 download per organism via NCBI Datasets (cached under
# outputs/genomes_by_organism/), matched locally against family_labels.csv by
# exact NCBI GeneID -> gene symbol -> family keyword, in that preference
# order. Far fewer network calls and far fewer chances of a fuzzy per-gene
# text search grabbing the wrong record than the old per-gene Entrez
# approach.

def ensure_datasets_cli():
    if shutil.which("datasets") is None:
        sys.exit(
            "[!] NCBI 'datasets' command-line tool not found on PATH -- needed to "
            "download genomes for promoter extraction.\n"
            "    Install it (precompiled binary, no build needed):\n"
            "        conda install -c conda-forge ncbi-datasets-cli\n"
            "    or see: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
        )


def extract_strain_tokens(organism):
    """Best-effort: strip the leading Genus + species (first two words) and
    treat whatever's left as strain/isolate identifiers to match against.
    "Plasmodium knowlesi strain H" -> ["H"]; "Toxoplasma gondii
    GAB2-2007-GAL-DOM2" -> ["GAB2-2007-GAL-DOM2"]. Two-word species-only
    queries (e.g. "Plasmodium vivax") correctly yield no tokens -- there's
    no strain to match, so the species reference IS the right answer."""
    parts = organism.strip().split()
    tail = parts[2:] if len(parts) > 2 else []
    tokens = [t for t in tail if t.lower() not in ("strain", "isolate", "str.", "str")]
    return tokens or tail


def resolve_strain_specific_accession(organism):
    """
    Lists EVERY assembly NCBI has for this taxon and picks the one that
    actually matches the requested strain, instead of trusting
    `--reference` to do it.

    WHY THIS EXISTS (bug found 2026-09): `datasets download genome taxon
    "<name>" --reference` silently resolves to whichever assembly NCBI has
    flagged as THE reference for the SPECIES that taxon belongs to -- and
    per NCBI's own docs, there is normally exactly ONE reference assembly
    per species, not one per strain. For a strain-level query like
    "Toxoplasma gondii ARI", there is usually no reference-flagged
    assembly under the ARI-specific taxon node, so the CLI's taxon match
    falls back up to the species and quietly hands back the SAME
    species-wide reference genome (e.g. ME49 for every T. gondii strain
    query) no matter which strain string was requested. Every
    differently-named strain folder then ends up holding an identical
    genome, so GFF GeneIDs for all but the one real reference strain never
    match family_labels.csv's ncbi_gene_id values (see
    build_gene_lookup/match_genes_in_gff) -- this is what was silently
    discarding ~86% of genes upstream of promoter extraction.

    Returns (accession, match_quality) where match_quality is one of:
      "strain_match"              -- found metadata containing the strain token(s)
      "only_assembly"             -- taxon has exactly one assembly; unambiguous
      "species_reference_fallback"-- multiple candidates, none strain-specific;
                                      used whichever is flagged reference/representative
      "arbitrary_fallback"        -- multiple candidates, none flagged either
      (None, None)                -- `datasets summary` returned nothing usable
    """
    strain_tokens = extract_strain_tokens(organism)

    cmd = ["datasets", "summary", "genome", "taxon", organism, "--as-json-lines"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        return None, None

    candidates = []
    for line in result.stdout.strip().splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        acc = rec.get("accession")
        if not acc:
            continue
        org = rec.get("organism", {}) or {}
        infra = org.get("infraspecific_names", {}) or {}
        searchable = " ".join(str(v) for v in infra.values() if v).lower()
        searchable += " " + str(org.get("organism_name", "")).lower()
        candidates.append((acc, searchable, rec))

    if not candidates:
        return None, None

    if strain_tokens:
        for acc, searchable, _rec in candidates:
            if any(tok.lower() in searchable for tok in strain_tokens):
                return acc, "strain_match"

    if len(candidates) == 1:
        return candidates[0][0], "only_assembly"

    for acc, _searchable, rec in candidates:
        category = ((rec.get("assembly_info", {}) or {}).get("refseq_category", "") or "")
        if "reference" in category.lower() or "representative" in category.lower():
            return acc, "species_reference_fallback"

    return candidates[0][0], "arbitrary_fallback"


def download_genome_for_organism(organism, dest_dir):
    """Downloads the correct genome assembly (FASTA + GFF3) for `organism`
    via NCBI Datasets, into dest_dir/genome.fna and dest_dir/genome.gff.
    Skips the download if both files already exist.

    Resolves to a SPECIFIC assembly accession first (see
    resolve_strain_specific_accession) rather than trusting a bare
    `taxon ... --reference` query, which silently collapses distinct
    strain-level queries onto the one species-wide reference assembly --
    see that function's docstring for the full bug writeup."""
    fna_path = os.path.join(dest_dir, "genome.fna")
    gff_path = os.path.join(dest_dir, "genome.gff")
    if os.path.exists(fna_path) and os.path.exists(gff_path):
        print(f"    ✓ Already downloaded -> {dest_dir}")
        return fna_path, gff_path

    os.makedirs(dest_dir, exist_ok=True)
    zip_path = os.path.join(dest_dir, "ncbi_dataset.zip")

    accession, match_quality = resolve_strain_specific_accession(organism)

    if accession:
        if match_quality == "strain_match":
            print(f"    Resolved '{organism}' -> strain-matched assembly {accession}")
        elif match_quality == "only_assembly":
            print(f"    Resolved '{organism}' -> only available assembly {accession}")
        else:
            print(f"    [!] No strain-specific assembly metadata found for '{organism}' -- "
                  f"falling back to the species-level reference/representative assembly "
                  f"{accession} ({match_quality}). GeneIDs for this specific strain may "
                  f"NOT match this genome's annotation -- check the match_method column "
                  f"in promoter_fetch_log.csv for this organism afterward.")
        # Record what was actually downloaded, right next to the genome
        # files, so a mismatch is inspectable later without re-running
        # anything.
        with open(os.path.join(dest_dir, "assembly_used.txt"), "w") as f:
            f.write(f"requested_organism={organism}\naccession={accession}\n"
                    f"match_quality={match_quality}\n")
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome,gff3",
            "--filename", zip_path
        ]
    else:
        print(f"    [!] 'datasets summary' returned nothing for '{organism}' -- falling back "
              f"to the old taxon+--reference query (may silently grab the wrong strain).")
        cmd = [
            "datasets", "download", "genome", "taxon", organism,
            "--reference",
            "--include", "genome,gff3",
            "--filename", zip_path
        ]

    print(f"    Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(zip_path):
        print(f"    [!] Genome download failed for '{organism}'.")
        if result.stderr:
            print(f"        {result.stderr.strip()[-500:]}")
        return None, None

    extract_dir = os.path.join(dest_dir, "_extracted")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)

    fna_hits = glob.glob(os.path.join(extract_dir, "**", "*.fna"), recursive=True)
    gff_hits = glob.glob(os.path.join(extract_dir, "**", "*.gff"), recursive=True)
    if not fna_hits or not gff_hits:
        print(f"    [!] Download for '{organism}' didn't contain both a genome FASTA and GFF3.")
        shutil.rmtree(extract_dir, ignore_errors=True)
        return None, None

    shutil.copy(fna_hits[0], fna_path)
    shutil.copy(gff_hits[0], gff_path)
    shutil.rmtree(extract_dir, ignore_errors=True)
    os.remove(zip_path)

    print(f"    ✓ Downloaded and cached -> {dest_dir}")
    return fna_path, gff_path


def parse_dbxref_geneid(attrs):
    """Extract the NCBI GeneID from a GFF3 Dbxref attribute like
    'Dbxref=GeneID:672,Genbank:NP_001.1' -> '672'. Returns None if absent."""
    dbxref = attrs.get("Dbxref", "")
    for token in dbxref.split(","):
        token = token.strip()
        if token.upper().startswith("GENEID:"):
            return token.split(":", 1)[1].strip()
    return None


def extract_promoter_seq_genomic(chr_seq, start, end, strand, upstream_bp, downstream_bp):
    """Strand-aware promoter extraction, TSS-relative (upstream_bp before the
    TSS through downstream_bp past it). start/end are 1-based inclusive GFF
    coordinates in reference (+) orientation."""
    chr_len = len(chr_seq)
    if strand == "+":
        tss = start - 1
        region_start = max(0, tss - upstream_bp)
        region_end = min(chr_len, tss + downstream_bp)
        return str(chr_seq[region_start:region_end])
    else:
        tss = end
        region_start = max(0, tss - downstream_bp)
        region_end = min(chr_len, tss + upstream_bp)
        return str(chr_seq[region_start:region_end].reverse_complement())


def build_gene_lookup(labels_df, organism):
    """For one organism, build lookup dicts from family_labels.csv rows:
    ncbi_gene_id -> [seq_id, ...], and gene symbol (lowercased) -> [seq_id, ...].

    CORRECTNESS FIX (found 2026-09): these are lists, not single values.
    family_labels.csv legitimately has multiple seq_id rows -- different
    UniProt/RefSeq protein accessions for the same underlying gene (e.g. an
    XP_ RefSeq record and an EPT GenBank record for one gene) -- sharing the
    SAME ncbi_gene_id. A plain dict here means the second row's seq_id
    silently overwrites the first's, permanently orphaning it even though
    its gene really was found in the GFF. This is exactly why
    Toxoplasma_gondii_ME49 (132 rows, only 66 distinct gene_ids) produced
    only 66 promoters even with a correctly-downloaded genome and complete
    GeneIDs: half its rows were shadowed out of the old single-value dict.
    A promoter belongs to the gene locus, not to a specific duplicate
    protein record, so every seq_id sharing a gene_id/symbol should get the
    same promoter -- see match_genes_in_gff below."""
    org_rows = labels_df[labels_df["organism"].astype(str) == organism]
    by_geneid, by_symbol = defaultdict(list), defaultdict(list)
    for _, row in org_rows.iterrows():
        gid = str(row.get("ncbi_gene_id", "")).strip()
        if gid and gid.lower() != "nan":
            by_geneid[gid].append(row["seq_id"])
        gene_sym = str(row.get("gene", "")).strip().lower()
        if gene_sym and gene_sym != "nan":
            by_symbol[gene_sym].append(row["seq_id"])
    return org_rows, by_geneid, by_symbol


def match_genes_in_gff(gff_path, fna_path, labels_df, organism, family_keyword,
                        upstream_bp, downstream_bp):
    """Scans one organism's downloaded GFF for gene features, matches each
    against family_labels.csv rows for that organism (GeneID -> symbol ->
    keyword, in that preference order), and extracts the promoter for every
    match. Returns a list of result dicts.

    NOTE: a single GFF gene feature can now satisfy MULTIPLE seq_ids at once
    (every family_labels.csv row sharing that gene's ncbi_gene_id or
    symbol) -- see build_gene_lookup's docstring for why that's correct,
    not a widened match."""
    org_rows, by_geneid, by_symbol = build_gene_lookup(labels_df, organism)
    wanted_seq_ids = set(org_rows["seq_id"])
    matched_seq_ids = set()

    print(f"    Loading genome sequence(s) from {fna_path} ...")
    genome_seqs = SeqIO.to_dict(SeqIO.parse(fna_path, "fasta"))

    results = []
    with open(gff_path, "r") as gff_handle:
        for line in gff_handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature_type, start_str, end_str, score, strand, phase, attributes_raw = parts
            if feature_type.lower() not in ("gene", "protein_coding_gene", "pseudogene"):
                continue

            attrs = parse_gff_attributes(attributes_raw)
            gff_gene_id = parse_dbxref_geneid(attrs)
            gff_symbol = (attrs.get("Name") or attrs.get("gene") or attrs.get("locus_tag") or "").strip().lower()
            description = (attrs.get("description") or attrs.get("Note") or "").lower()

            candidate_seq_ids, match_method = [], None
            if gff_gene_id and gff_gene_id in by_geneid:
                candidate_seq_ids = by_geneid[gff_gene_id]
                match_method = "exact_geneid"
            elif gff_symbol and gff_symbol in by_symbol:
                candidate_seq_ids = by_symbol[gff_symbol]
                match_method = "gene_symbol"
            elif family_keyword and family_keyword.lower() in description:
                unmatched = wanted_seq_ids - matched_seq_ids
                if unmatched:
                    # Keyword fallback stays single-assignment: it's already
                    # an arbitrary/ambiguous pick, so widening it to "every
                    # remaining unmatched seq_id" would make it worse, not
                    # better.
                    candidate_seq_ids = [next(iter(unmatched))]
                    match_method = "family_keyword_fallback"

            new_seq_ids = [sid for sid in candidate_seq_ids
                           if sid in wanted_seq_ids and sid not in matched_seq_ids]
            if not new_seq_ids:
                continue
            if seqid not in genome_seqs:
                continue

            start, end = int(start_str), int(end_str)
            chr_seq = genome_seqs[seqid].seq
            promoter_seq = extract_promoter_seq_genomic(chr_seq, start, end, strand, upstream_bp, downstream_bp)
            if not promoter_seq:
                continue

            for seq_id in new_seq_ids:
                matched_seq_ids.add(seq_id)
                results.append({
                    "seq_id": seq_id, "gene_id": gff_gene_id or "", "status": "ok",
                    "sequence": promoter_seq, "organism": organism, "match_method": match_method
                })

    for seq_id in (wanted_seq_ids - matched_seq_ids):
        results.append({"seq_id": seq_id, "status": "gene_not_found_in_gff",
                         "organism": organism, "match_method": ""})
    return results


def fetch_and_write_promoters_via_genome(labels_path="outputs/family_labels.csv",
                                          upstream_bp=PROMOTER_UPSTREAM_BP,
                                          downstream_bp=PROMOTER_DOWNSTREAM_BP):
    """Promoter path for Mode 2: reads family_labels.csv, downloads one
    genome+GFF per distinct organism, matches genes locally, and writes the
    same promoter output namespace (outputs/promoter_sequences.fasta,
    outputs/promoter_fetch_log.csv, outputs/promoters_by_organism/<Org>/)
    that Mode 1 also produces."""
    import pandas as pd

    if not os.path.exists(labels_path):
        print(f"[!] {labels_path} not found -- skipping promoter fetch.")
        return

    ensure_datasets_cli()

    labels_df = pd.read_csv(labels_path)
    if "ncbi_gene_id" in labels_df.columns:
        # Avoid pandas silently upcasting a partially-blank int column to
        # float64 (1001 -> "1001.0"), which would break exact GeneID matches
        # against the GFF's plain integer strings.
        labels_df["ncbi_gene_id"] = (
            labels_df["ncbi_gene_id"].astype(str)
            .str.replace(r"\.0$", "", regex=True)
            .replace("nan", "")
        )
    if "organism" not in labels_df.columns:
        print(f"[!] {labels_path} has no 'organism' column -- skipping promoter fetch.")
        return

    family_keyword = str(labels_df["family"].iloc[0]) if "family" in labels_df.columns and len(labels_df) else ""
    organisms = sorted(labels_df["organism"].dropna().astype(str).unique())

    print(f"\n[*] Fetching promoters for {len(labels_df)} sequence(s) across "
          f"{len(organisms)} organism(s), via genome download + local GFF matching.")
    print(f"    Region: {upstream_bp}bp upstream to {downstream_bp}bp downstream of TSS\n")

    all_results = []
    os.makedirs(GENOME_ROOT, exist_ok=True)

    for organism in organisms:
        print(f"{'=' * 60}\nOrganism: {organism}\n{'=' * 60}")
        org_dir = os.path.join(GENOME_ROOT, sanitize_organism_name(organism))

        fna_path, gff_path = download_genome_for_organism(organism, org_dir)
        if not fna_path or not gff_path:
            for seq_id in labels_df.loc[labels_df["organism"].astype(str) == organism, "seq_id"]:
                all_results.append({"seq_id": seq_id, "status": "genome_download_failed",
                                     "organism": organism, "match_method": ""})
            continue

        org_results = match_genes_in_gff(gff_path, fna_path, labels_df, organism,
                                          family_keyword, upstream_bp, downstream_bp)
        n_ok = sum(1 for r in org_results if r["status"] == "ok")
        print(f"    ✓ Matched and extracted {n_ok}/{len(org_results)} gene(s) for this organism.")
        all_results.extend(org_results)

    _write_promoter_outputs(all_results)


def split_local_promoters_by_organism(labels_path="outputs/family_labels.csv",
                                       promoter_fasta_path="outputs/promoter_sequences.fasta"):
    """Promoter path for Mode 1: promoters were already extracted directly
    from the local GFF (fast, no download needed, and already correctly
    matched by construction). This just re-reads what write_outputs() wrote
    and splits/logs it into the same output namespace the genome-download
    path uses, so downstream tooling doesn't need to special-case the mode
    that produced outputs/promoter_sequences.fasta."""
    import pandas as pd

    if not (os.path.exists(labels_path) and os.path.exists(promoter_fasta_path)):
        return

    labels_df = pd.read_csv(labels_path)
    seq_to_org = dict(zip(labels_df["seq_id"], labels_df["organism"].astype(str)))

    results = []
    for rec in SeqIO.parse(promoter_fasta_path, "fasta"):
        org = seq_to_org.get(rec.id)
        results.append({
            "seq_id": rec.id, "gene_id": "", "status": "ok" if org else "organism_unmatched",
            "sequence": str(rec.seq), "organism": org or "", "match_method": "local_gff_direct"
        })

    _write_promoter_outputs(results, write_master_fasta=False)


def _write_promoter_outputs(results, write_master_fasta=True):
    """Shared writer for both promoter-fetch paths: the combined
    outputs/promoter_fetch_log.csv (uniform columns regardless of mode) and
    the outputs/promoters_by_organism/<Organism>/ per-organism split.
    Mode 1 already wrote outputs/promoter_sequences.fasta itself via
    write_outputs(), so write_master_fasta=False there to avoid clobbering it
    with a possibly-reordered rewrite."""
    os.makedirs("outputs", exist_ok=True)
    n_ok = sum(1 for r in results if r["status"] == "ok")

    if write_master_fasta:
        with open("outputs/promoter_sequences.fasta", "w") as f:
            for r in results:
                if r["status"] == "ok":
                    f.write(f">{r['seq_id']}\n{r['sequence']}\n")

    with open("outputs/promoter_fetch_log.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["seq_id", "gene_id", "status", "organism", "match_method"])
        for r in results:
            writer.writerow([r["seq_id"], r.get("gene_id", ""), r["status"],
                              r.get("organism", ""), r.get("match_method", "")])

    print(f"\n[✓] Promoters: {n_ok}/{len(results)} sequence(s) resolved.")
    print("    FASTA -> outputs/promoter_sequences.fasta")
    print("    Log   -> outputs/promoter_fetch_log.csv")

    org_root = os.path.join("outputs", "promoters_by_organism")
    os.makedirs(org_root, exist_ok=True)
    by_organism = {}
    for r in results:
        if r["status"] != "ok":
            continue
        org_key = sanitize_organism_name(r.get("organism", ""))
        by_organism.setdefault(org_key, []).append(r)

    for org_key, org_results in by_organism.items():
        organism_dir = os.path.join(org_root, org_key)
        os.makedirs(organism_dir, exist_ok=True)
        with open(os.path.join(organism_dir, "promoter_sequences.fasta"), "w") as f:
            for r in org_results:
                f.write(f">{r['seq_id']}\n{r['sequence']}\n")
        with open(os.path.join(organism_dir, "promoter_fetch_log.csv"), "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["seq_id", "gene_id", "status", "organism", "match_method"])
            for r in org_results:
                writer.writerow([r["seq_id"], r.get("gene_id", ""), r["status"],
                                  r.get("organism", ""), r.get("match_method", "")])

    print(f"    Organism-wise outputs -> {org_root}/<Organism>/")
    for org_key, org_results in sorted(by_organism.items()):
        print(f"      {org_key}: {len(org_results)} promoter(s)")


# ==============================================================================
# OUTPUT WRITER
# ==============================================================================

# ==============================================================================
# PROTEIN vs. NUCLEOTIDE DETECTION + DATA MANIFEST
#
# Local mode (process_local_genome_and_gff) writes the RAW GENOMIC sequence
# straight from the GFF coordinates into "Sequence" -- there is no
# .translate() call anywhere in this file, so despite the module docstring
# calling family_sequences.fasta a "master protein FASTA", local mode's
# output is actually nucleotide. NCBI mode (fetch_ncbi_protein_online) pulls
# from the "protein" Entrez database directly, so its output is genuinely
# protein.
#
# Rather than hardcoding "local mode == no protein" (which breaks silently
# the moment this script is extended to translate CDS locally, or to accept
# a pre-translated local protein FASTA), this checks the ACTUAL sequence
# content written to fasta_path, so downstream stages get a manifest that's
# always correct regardless of which code path produced the data.
# ==============================================================================

NUCLEOTIDE_ALPHABET = set("ACGTUN")
DATA_MANIFEST_PATH = "outputs/.pipeline_data_manifest.json"


def detect_sequence_type(fasta_path, sample_size=200):
    """
    Returns ("protein", fraction_non_nucleotide) or ("nucleotide", fraction)
    or ("empty", 0.0) by sampling up to sample_size sequences from
    fasta_path and checking what fraction of residues fall outside
    A/C/G/T/U/N. Real protein sequences are overwhelmingly NOT in that set
    (most residues are e.g. L/E/S/K/etc.); nucleotide sequences are almost
    entirely within it.
    """
    if not os.path.exists(fasta_path):
        return "empty", 0.0

    total_chars = 0
    non_nucleotide_chars = 0
    n_seqs_checked = 0

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                if line.startswith(">"):
                    n_seqs_checked += 1
                    if n_seqs_checked > sample_size:
                        break
                continue
            total_chars += len(line)
            non_nucleotide_chars += sum(1 for c in line.upper() if c not in NUCLEOTIDE_ALPHABET)

    if total_chars == 0:
        return "empty", 0.0

    fraction_non_nt = non_nucleotide_chars / total_chars
    # A handful of ambiguity codes (R/Y/etc.) can sneak into real nucleotide
    # data, so use a permissive threshold rather than requiring 100% purity.
    seq_type = "nucleotide" if fraction_non_nt < 0.05 else "protein"
    return seq_type, fraction_non_nt


def write_data_manifest(fasta_path, mode, manifest_path=DATA_MANIFEST_PATH,
                         known_protein_counts=None):
    """
    Writes outputs/.pipeline_data_manifest.json so downstream pipeline
    stages (master_script.py / master_script_2.py) know, without guessing,
    whether real protein data is available for this run -- see those
    scripts' load_data_manifest() / has_real_protein handling.

    known_protein_counts, when given, is (n_real_protein, n_total) from the
    fetcher's OWN per-gene bookkeeping (currently: local mode's protein-FASTA
    ID matching -- see process_local_genome_and_gff's _is_real_protein
    flag). When available this is authoritative and takes priority over the
    content-sniffing heuristic below, because it's a certainty rather than
    an inference -- and critically, it catches a MIXED dataset (some genes
    matched to a real protein, others fell back to nucleotide) that a
    single aggregate content check across the whole FASTA could misjudge
    either way, depending on which type happens to dominate by sequence
    count.
    """
    seq_type, fraction_non_nt = detect_sequence_type(fasta_path)
    has_real_protein = (seq_type == "protein")
    is_mixed = False

    if known_protein_counts is not None:
        n_real_protein, n_total = known_protein_counts
        if n_total > 0:
            has_real_protein = n_real_protein > 0
            is_mixed = 0 < n_real_protein < n_total
            seq_type = "protein" if n_real_protein == n_total else (
                "mixed" if is_mixed else "nucleotide"
            )

    manifest = {
        "fasta_path": fasta_path,
        "fetch_mode": mode,                       # "local" or "ncbi"
        "detected_sequence_type": seq_type,        # "protein" / "nucleotide" / "mixed" / "empty"
        "fraction_non_nucleotide_chars": round(fraction_non_nt, 4),
        "has_real_protein": has_real_protein,
        "written_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    if known_protein_counts is not None:
        manifest["known_protein_gene_count"] = known_protein_counts[0]
        manifest["known_total_gene_count"] = known_protein_counts[1]

    os.makedirs(os.path.dirname(manifest_path) or ".", exist_ok=True)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\n[*] Data manifest -> {manifest_path}")
    print(f"    Detected sequence type : {seq_type}")
    print(f"    Real protein data      : {'YES' if has_real_protein else 'NO'}")

    if is_mixed:
        n_real_protein, n_total = known_protein_counts
        print(
            f"\n[!] MIXED dataset: {n_real_protein}/{n_total} gene(s) matched a real "
            f"protein sequence; the other {n_total - n_real_protein} fall back to their "
            f"nucleotide gene sequence in the SAME {fasta_path}."
            "\n    has_real_protein is set to YES so downstream protein-analysis stages "
            "still run, but they will run on the WHOLE file -- the nucleotide-fallback "
            "entries will be silently scored as if they were protein too (this is a "
            "known limitation, not caught automatically: there is currently no "
            "per-sequence protein/nucleotide tag carried into family_labels.csv for "
            "those stages to filter on). If that fraction is large, treat protein-stage "
            "results with caution, or exclude organisms with poor protein-FASTA "
            "coverage from this run."
        )

    if not has_real_protein:
        if mode == "local":
            reason = ("no companion protein FASTA was found (or none of its IDs matched "
                       "the GFF's gene IDs) for these genes, so nucleotide gene sequences "
                       "were used instead -- see find_matching_protein_fasta() if you do "
                       "have local protein FASTAs and expected them to be picked up")
        else:
            reason = "unexpected for NCBI protein mode -- double-check the fetch above"
        print(
            f"\n[!] No real protein sequences detected for this run ({reason})."
            "\n    Protein-only analysis stages (protein property profiling, "
            "domain architecture) will need to be skipped."
        )
        print(
            "    To continue the automated pipeline past that, on BOTH the "
            "combined and genus/organism-wise master scripts, run:"
        )
        print("        python master_script.py --from 5")
        print("        python master_script_2.py --from 5")
        print(
            "    (both scripts read this manifest automatically and will "
            "skip protein-only stages and enforce Stage 05+ on their own -- "
            "see NO_PROTEIN_MIN_STAGE near the top of each. Adjust that "
            "constant if your protein-dependent stages sit at a different "
            "stage number.)"
        )

    return manifest


def write_outputs(records, family_name, fasta_path, labels_path, architecture_path, promoter_records=None, fetch_mode="local"):
    seen = set()
    arch_rows = []
    n_written = 0

    org_dir = os.path.join(os.path.dirname(fasta_path), "by_organism")
    os.makedirs(org_dir, exist_ok=True)
    
    org_files = {}

    with open(fasta_path, "w") as fasta, open(labels_path, "w", newline="") as labcsv:
        writer = csv.writer(labcsv)
        writer.writerow(["seq_id", "family", "length", "organism", "gene",
                         "domain_trimmed", "ensembl_transcript_id", "ncbi_gene_id"])
        
        for row in records:
            acc = row["Entry"]
            if acc in seen:
                continue
            seen.add(acc)
            
            seq = row["Sequence"].strip()
            if not seq:
                continue
                
            org = row.get("Organism", "unknown")
            gene_raw = row.get("Gene Names", "")
            gene = gene_raw.split()[0] if gene_raw else acc
            
            transcript_id = row.get("_ensembl_transcript") or ""
            seq_id = make_seq_id(acc, gene)
            
            fasta_entry = f">{seq_id}\n{seq}\n"
            
            # Write to Master FASTA
            fasta.write(fasta_entry)
            
            # Write to Individual Organism FASTA
            clean_org = sanitize_organism_name(org)
            if clean_org not in org_files:
                org_filepath = os.path.join(org_dir, f"{clean_org}.fasta")
                org_files[clean_org] = open(org_filepath, "w")
            
            org_files[clean_org].write(fasta_entry)

            # Write Metadata CSV
            writer.writerow([seq_id, family_name, len(seq), org, gene,
                             row["_domain_trimmed"], transcript_id, row.get("_ncbi_gene_id", "")])
            n_written += 1

            full_length = row["_full_length"]
            for dom_name, dom_start, dom_end in row["_all_domains"]:
                arch_rows.append([seq_id, full_length, dom_name, dom_start, dom_end])

    # Close all open organism file handles
    for f in org_files.values():
        f.close()

    # Write Domain Architecture CSV
    with open(architecture_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["seq_id", "protein_length", "domain_name", "start", "end"])
        writer.writerows(arch_rows)

    # Write Deduplicated Promoters
    if promoter_records:
        promoter_fasta_path = "outputs/promoter_sequences.fasta"
        seen_promoters = set()
        unique_count = 0
        
        with open(promoter_fasta_path, "w") as p_fasta:
            for p_id, p_seq in promoter_records:
                if p_id not in seen_promoters:
                    seen_promoters.add(p_id)
                    p_fasta.write(f">{p_id}\n{p_seq}\n")
                    unique_count += 1
                    
        print(f"  • Promoter FASTA -> {promoter_fasta_path} ({unique_count} unique gene promoters written)")

    print(f"\n[✓] Outputs generated successfully:")
    print(f"  • Master FASTA     → {fasta_path}")
    print(f"  • Organism FASTAs → {org_dir}/ ({len(org_files)} organism file(s) created)")
    print(f"  • Labels          → {labels_path}")
    print(f"  • Architecture    → {architecture_path}")

    # Only local mode's records carry the _is_real_protein bookkeeping (set
    # by process_local_genome_and_gff when a companion protein FASTA was
    # discovered and ID-matched). NCBI mode's records never set this key at
    # all -- its content is unconditionally real protein already, and the
    # content-sniffing heuristic in write_data_manifest handles that fine
    # on its own, so we deliberately do NOT pass known_protein_counts for
    # fetch_mode == "ncbi" (that would wrongly report 0/N there, since
    # NCBI records simply never carry this key to count as True).
    known_protein_counts = None
    if fetch_mode == "local":
        rows_with_flag = [r for r in records if "_is_real_protein" in r]
        if rows_with_flag:
            n_real_protein = sum(1 for r in rows_with_flag if r["_is_real_protein"])
            known_protein_counts = (n_real_protein, len(rows_with_flag))

    write_data_manifest(fasta_path, mode=fetch_mode, known_protein_counts=known_protein_counts)

    return n_written


# ==============================================================================
# MAIN ENTRY POINT
# ==============================================================================

if __name__ == "__main__":
    os.makedirs("outputs", exist_ok=True)
    
    print("=" * 65)
    print("     GENE FAMILY & PROMOTER DATA EXTRACTION TOOL       ")
    print("     (fetches proteins AND promoters in one run)               ")
    print("=" * 65)
    print(" Select Data Retrieval Mode:")
    print("  [1] Extract from Local Genome FASTA & GFF Files (Default)")
    print("      -> promoters extracted directly from the same local GFF pass")
    print("  [2] Fetch from NCBI Database Online")
    print("      -> promoters via NCBI genome download + local GFF match")
    print("=" * 65)

    choice = input("Enter choice [1 or 2, Default: 1]: ").strip()

    if choice == "2":
        print("\n--- ONLINE MODE: NCBI PROTEIN FETCH ---")
        org_list = get_organism_list_from_user(default_orgs=["Plasmodium", "Toxoplasma", "Babesia"])

        family_input = input("\nEnter Gene Family / Protein Name Keyword [default: AP2 domain]: ").strip()
        if not family_input:
            family_input = "AP2 domain"

        records, family_label = fetch_ncbi_protein_online(family_input, org_list)

        if records:
            write_outputs(records,
                          family_label,
                          "outputs/family_sequences.fasta",
                          "outputs/family_labels.csv",
                          "outputs/domain_architecture.csv",
                          fetch_mode="ncbi")
            # No promoter_records from this mode -- fetch them separately by
            # downloading each organism's genome+GFF and matching locally.
            # Wrapped so a promoter-fetch failure (e.g. `datasets` CLI
            # missing, or a genome download failing for one organism)
            # doesn't discard the protein fetch that already succeeded.
            try:
                fetch_and_write_promoters_via_genome()
            except SystemExit as e:
                print(f"\n[!] Promoter fetch skipped: {e}")
            except Exception as e:
                print(f"\n[!] Promoter fetch failed unexpectedly: {e}")

    else:
        print("\n--- LOCAL MODE: GENOME & GFF EXTRACTION ---")
        print("Every FASTA/GFF pair found under each directory you give (subfolders "
              "are searched too) is treated as its own organism/strain -- organism "
              "names are read from each genome FASTA's own header when present, "
              "not guessed from folder names.")

        root_dirs = get_local_root_dirs_from_user(default_dir="Genomes")

        pairs = []
        seen_pair_keys = set()
        print()
        for root_dir in root_dirs:
            if not os.path.isdir(root_dir):
                print(f"  [!] Directory '{root_dir}' does not exist -- skipping.")
                continue

            found = discover_gff_fasta_pairs(root_dir)

            # Merge across directories without double-counting a pair that
            # happens to be reachable from more than one provided path
            # (e.g. one path is a subfolder of another).
            new_found = []
            for p in found:
                key = (os.path.abspath(p["gff"]), os.path.abspath(p["fasta"]))
                if key in seen_pair_keys:
                    continue
                seen_pair_keys.add(key)
                new_found.append(p)

            n_dupes = len(found) - len(new_found)
            print(f"  [*] {root_dir}: found {len(new_found)} genome FASTA + GFF pair(s)"
                  + (f" ({n_dupes} already seen from another directory, skipped)" if n_dupes else "")
                  + ".")
            pairs.extend(new_found)

        if not pairs:
            sys.exit(f"[!] No genome FASTA + GFF pairs found under any of: {', '.join(root_dirs)}.")

        chosen_pairs = get_pairs_from_user(pairs)
        if not chosen_pairs:
            sys.exit("[!] No genome/GFF pairs selected.")

        gene_keyword = input("\nEnter Gene Family / Domain Keyword to filter (e.g., Actin, AP2, Kinase): ").strip()
        while not gene_keyword:
            gene_keyword = input("Keyword cannot be empty. Please enter a target keyword: ").strip()

        # Locked: every downstream stage (06a/06b/06c/07i/05g/010/...) assumes
        # PROMOTER_LEN bp upstream-only promoters, so this is no longer a prompt.
        prom_len = PROMOTER_LEN
        print(f"[i] Promoter length locked to {prom_len} bp upstream of the TSS (pipeline_config.PROMOTER_LEN).")

        records, promoter_records, family_label = process_local_genomes_multi_org(
            chosen_pairs, gene_keyword, prom_len)

        if records:
            write_outputs(records,
                          family_label,
                          "outputs/family_sequences.fasta",
                          "outputs/family_labels.csv",
                          "outputs/domain_architecture.csv",
                          promoter_records=promoter_records,
                          fetch_mode="local")
            # Promoters were already extracted directly from the local GFF
            # above (correctly keyed via make_seq_id -- see
            # process_local_genome_and_gff). Split what write_outputs() just
            # wrote into the same per-organism namespace the online mode uses.
            split_local_promoters_by_organism()
        else:
            print("[!] No matching features found in any organism's GFF annotations.")

    print("\nNext step: Run '02_align_with_mafft_for_all.py'")