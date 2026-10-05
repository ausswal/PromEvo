# Comparative Promoter/Gene-Family Analysis Pipeline

A comparative genomics pipeline for analyzing one gene family's promoters
and coding sequences across multiple species/genera -- phylogenetics
(gene-level, genus-level, combined, and an OrthoFinder-based reference
species tree), motif discovery (MEME/TOMTOM), core promoter element
analysis, CpG islands, DNA structural properties, and comparative HTML
dashboards.

## Environment setup

Two equivalent ways to get the exact same environment:

**Option A -- one command (recommended for first-time setup):**
```bash
bash run_pipeline.sh
```
This creates the `gene_promoter_env` conda environment, downloads the
shared JASPAR/UniProbe motif reference databases into `databases/`
(one-time, reused by every future run), then launches Stage 01's
interactive prompt. If you're not ready to start a real run yet, it's
safe to `Ctrl+C` once it reaches that prompt -- the environment and
databases are already fully set up by that point.

**Option B -- from the environment spec directly:**
```bash
conda env create -f environment.yml
conda activate gene_promoter_env
python 00_setup_motif_databases.py
```
Use this if you want the environment without immediately launching the
pipeline, or if you're recreating the environment on a new machine and
want the standard `conda env create` workflow.

## Verifying the setup

Before starting a real run (which can take hours), check that packages,
CLI tools, and every script's syntax are all in order:
```bash
conda activate gene_promoter_env
python environment_check.py
```
This never touches `outputs/` or runs any pipeline stage -- it only
checks that a real run would have what it needs.

## Running the pipeline

```bash
conda activate gene_promoter_env
python master_script_8.py
```

Useful flags (see `python master_script_8.py --list` for the full stage
list with completion status):
- `--resume` -- continue automatically after the last successfully
  completed stage.
- `--from N` -- jump straight to stage `N`, skipping everything before it.
- `--reset` -- clear saved progress and start fresh next time.

## Local vs. NCBI-fetched data

Stage 01 supports two modes:
- **NCBI mode** -- fetches sequences and promoters directly from NCBI for
  named organisms.
- **Local mode** -- scans one or more local directories of genome
  FASTA + GFF (optionally with a companion protein FASTA, e.g. VEuPathDB's
  `..._AnnotatedProteins.fasta`, for real protein-based analysis instead
  of nucleotide).

Whichever mode is used, the fetch stage writes
`outputs/.pipeline_data_manifest.json` recording whether real protein
data is available. Downstream stages read this automatically -- if only
nucleotide data is available, protein-only stages (protein property
profiling, the OrthoFinder reference species tree) are skipped cleanly
rather than producing meaningless output.

## Repository layout

- `master_script_8.py` -- orchestrates all pipeline stages in order.
- `run_pipeline.sh` -- one-command environment setup + launch.
- `environment.yml` -- conda environment specification (see above).
- `environment_check.py` -- pre-flight verification (packages/tools/syntax).
- `00_setup_motif_databases.py` -- one-time JASPAR/UniProbe database setup.
- `01e_new7_...py` through `010_...py` -- numbered pipeline stages.
- `011*_generate_comparative_dashboard.py` -- comparative HTML dashboard
  generators (genus-wise and species-wise, with an embedded on-the-fly
  download for the species-wise view).
- `Scripts_not_required/` -- superseded/experimental script versions kept
  for reference, not part of the active pipeline (excluded from version
  control -- see `.gitignore`).
- `outputs/` -- all generated results (excluded from version control;
  fully regenerable by re-running the pipeline).
- `databases/` -- shared JASPAR/UniProbe reference databases, set up once
  (excluded from version control; regenerable via
  `00_setup_motif_databases.py`).

## Notes on duplicate-looking script names

Some stages have more than one similarly-named script in this directory
(e.g. two `04b_...` variants, two `03c_...` variants). Only one of each
pair is referenced by `master_script_8.py`'s `STAGES` list at any given
time -- check that list directly if you're unsure which is the active
one before editing.
