# PromEvo

**Comparative gene-family and promoter analysis pipeline**

PromEvo analyses one gene family's promoters and coding sequences across
multiple species and genera: phylogenetics (gene-level, genus-level,
combined, and an OrthoFinder-based reference species tree), transcription
start site prediction, nucleotide composition, core promoter elements and
their synergism, motif discovery (MEME/TOMTOM), motif conservation and
distribution, CpG islands, DNA structural properties, non-B DNA motifs, and
comparative HTML dashboards.

## Installation

Once PromEvo is available on Bioconda:

```bash
conda create -n promevo -c conda-forge -c bioconda promevo
conda activate promevo
promevo-check
```

`promevo-check` verifies the Python packages, command-line tools (MEME Suite,
MAFFT, OrthoFinder, DIAMOND, MCL, FastTree, NCBI Datasets) and installed
pipeline files. It does not run any analysis or touch `outputs/`.

## Running the pipeline

Run PromEvo from a folder that contains your input genome folders. Results are
written to `./outputs`.

```bash
conda activate promevo
cd /path/to/your/project
promevo
```

On start-up PromEvo checks the environment, makes sure the motif reference
databases are present (see below), then launches Stage 01, which is
interactive. Run it in a normal terminal, not in the background.

Useful options (also available through `promevo-master`):

| Option | Effect |
|---|---|
| `--list` | List all stages with completion status, then exit |
| `--resume` | Continue after the last successfully completed stage |
| `--from N` | Start at stage `N`, skipping everything before it (Stage 01 is assumed done when `N > 1`) |
| `--reset` | Clear saved progress, then exit |

```bash
promevo-master --list
promevo-master --from 15       # rerun DNA structural properties onward
promevo-master --resume
```

To stop a run press **Ctrl+C**. Do not use Ctrl+Z, which only suspends it.

## Commands

| Command | Purpose |
|---|---|
| `promevo` | Full launcher: environment check, motif databases, then the pipeline |
| `promevo-master` | Run the master script directly (no database set-up step) |
| `promevo-check` | Verify the environment |
| `promevo-authors` / `promevo --authors` | Show author information |

## Pipeline stages

| Stage | Description |
|---|---|
| 01 | Retrieve gene family sequences and promoters (interactive) |
| 02 | Protein physicochemical properties (organism/genus-wise) |
| 03 | OrthoFinder reference species tree and phylogenetic trees |
| 04 | Multiple sequence alignment visualisation and conservation |
| 05 | Transcription start site prediction (upstream only) |
| 06 | Nucleotide composition analysis |
| 07 | Core promoter element analysis |
| 08 | Core promoter element synergism analysis |
| 09 | Motif discovery (MEME) and TOMTOM matching |
| 10 | Motif logos |
| 11 | Motif conservation analysis |
| 12 | Motif distribution (organism and genus) |
| 13 | CpG island finder |
| 14 | Gene-wise CpG plots |
| 15 | DNA structural properties |
| 16 | Non-B DNA motif analysis |
| 17 | Comparative HTML dashboard |

Each stage writes its log to `outputs/logs/`.

## Local vs. NCBI-fetched data

Stage 01 supports two modes:

- **Local mode (default)**: scans local directories of genome FASTA and GFF
  files.
- **NCBI mode**: fetches sequences and promoters directly from NCBI for named
  organisms (requires internet access).

Either way, Stage 01 writes `outputs/.pipeline_data_manifest.json`, recording
whether real protein data is available. Later stages read it automatically: if
only nucleotide data exists, protein-only stages (protein property profiling
and the OrthoFinder reference species tree) are skipped cleanly instead of
producing meaningless output.

## Motif reference databases

The JASPAR and UniProbe motif databases are not shipped with the package.
`promevo` downloads them on first run (internet access needed once) and reuses
them afterwards. They are stored in:

```
~/.local/share/comparative-gene-promoter/databases
```

To use another location, set `COMPARATIVE_PIPELINE_DATABASES` to the folder
before running, or set `XDG_DATA_HOME`.

## Repository layout

- `scripts/` : master script, numbered pipeline stages, launcher
  (`run_pipeline.sh`), environment checker (`environment_check.py`), motif
  database set-up (`00_setup_motif_databases.py`) and the dashboard generator
  (`011*_generate_comparative_dashboard.py`).
- `recipes/promevo/` : conda recipe (`meta.yaml`, `build.sh`).
- `AUTHORS.md`, `LICENSE.txt` : authorship and licence (MIT).
- `outputs/` : generated results, created when you run the pipeline in your
  project folder (not tracked in version control).

## Authors

Swarup Das¹, Subarna Thakur¹

¹ Department of Bioinformatics, University of North Bengal, Bagdogra,
Bairatisal, Darjeeling, West Bengal 734013, India

## License

MIT. See `LICENSE.txt`.
