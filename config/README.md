# Configuration

`config/config.yaml` has two parts.

| Key | Meaning |
|---|---|
| `scripts_dir` | Folder with the PromEvo stage scripts. A relative path is taken from the repository root (default `scripts`). Override on the command line with `--config scripts_dir=/absolute/path`. |
| `stages` | The 17 stages in order. Each entry has `id`, `script`, `after` (what must finish first; `databases` is the motif database setup), and optionally `guard: protein` (skipped automatically when Stage 01 found no real protein data) and `threads`. |

Stage 01 is interactive: it asks for the data-retrieval mode and inputs, then writes
`outputs/.pipeline_data_manifest.json`. Every later stage requires that file.

Finished steps are recorded as marker files in `outputs/.done/`. Delete `outputs/.done/1.done`
and `outputs/.done/{2..17}.done` to redo Stage 01 and everything after it.
