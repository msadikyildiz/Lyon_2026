# Genomic inputs and calculations

The repository includes all 149 per-sample mutation-call TSVs referenced by `reference/metadata_complete.csv`, plus seven ATEC and ten PLAC processed CSV exports. Their source archive is `Genomic_input_files.zip`. `input_manifest.json` records the archive checksum and each file's size, row count and SHA-256 checksum. Supplied files are preserved unchanged.

## Reproduce the tables

From the repository root, using the environment described in the main README:

```sh
python scripts/reproduce_genomics.py
```

This also runs as part of `python reproduce.py`, before the Source Data workbook is built. The script executes the original notebooks' data-processing cells in temporary directories, so their CSV exports cannot overwrite the supplied files. Notebook-specific filtering, allele handling, corrections in the executed data-processing cells and two-decimal frequency rounding are preserved. Display-label changes confined to plotting cells are not applied to these tables; those mappings remain in the original plotting code. Raw TSVs retain the original frequency precision. The code reproduces downstream processing of existing mutation calls; it does not align reads or call mutations.

| Output in `generated_tables/` | Rows | Figure coverage |
|---|---:|---|
| `MG_AMI` | 43 | Figure 2; Supplementary Figure 5 |
| `MG_LEV` | 64 | Figure 2; Supplementary Figure 5 |
| `MG_CEF` | 80 | Figure 2; Supplementary Figure 5 |
| `MG_CEF_R` | 56 | Supplementary Figure 3 |
| `PbEc_CEF` | 83 | Supplementary Figure 4 |
| `MG_untreated` | 38 | Supplementary Figure 6 |
| `PLAC` | 601 | Figure 4; Supplementary Figure 8 |
| `combined_lineages` | 963 | Supplementary Figure 9 |

Each table is exported as CSV and typed JSON. `generated_tables/index.json` maps all 12 figure notebooks to their output tables and records executed cell numbers and notebook hashes. The 17 supplied processed exports match the rerun calculations before the notebooks' later annotation/allele corrections. Every available saved endpoint row also matches after accounting for truncated display text, empty-field formatting and numerical serialization. `generated_tables/validation.json` records these checks.

The original input layout is retained:

```text
data/genomics/reference/metadata_complete.csv
data/genomics/out/<FolderDate>/<source_file>/output/output.gd.tsv
data/genomics/data/processed/traced_alleles/<lineage>/<population>.csv
```

`recovered_tables/` contains notebook-display exports used as validation references. Source Data use `generated_tables/`. The genomic plotting notebooks supply the selections and aliases used by scripts/reproduce_genomic_plots.py. All genomic numerical panels are rendered automatically. Plot matrices and coordinate/allele crosswalks are in plot_tables/. Distinct mutation sites remain separate even when they share a short label; variant identities are documented in `docs/DATA_CORRECTIONS.md`.

## References and run records

`Breseq_reference_and_run_records.zip` contains `data/reference.fasta`, `data/summary.json`, `output/summary.json` and `output/log.txt` for all 149 study runs. The original log and both summary JSON files remain under each dated run/sample directory. Two byte-distinct FASTAs are stored once under `reference/fasta/`; `run_records/manifest.json` maps every original archive path to its retained file and SHA-256. `run_records/runs.csv` maps each run to its command, options and archived GenBank reference.

All 141 MG1655 runs used `sequence-4.gb` (U00096.3); eight ATEC runs used `all_ATEC_annotated_contigs_w_locustag_w_genetags.gbk`. Both ATEC GenBanks are retained. For every run, the archived GenBank sequence matches the saved FASTA sequence, and the sequence length matches the summary JSON. Original `/work` reference paths no longer exist; sequence agreement does not establish byte identity of the original annotation files. Every recorded command supplied one R1-labelled FASTQ to breseq in polymorphism-prediction mode with 16 processors. This documents the analysis inputs, separately from the paired-end sequencing described in Methods. Saved options are retained in full, including the caller's 0.05 polymorphism-frequency cutoff; downstream notebook frequency filters remain separate.

Run `python scripts/verify_breseq_records.py` to verify every archive member, study-run identity and sequence match. This check also runs at the start of `reproduce.py`. The import option `--import-archive /path/to/Breseq_reference_and_run_records.zip` recreates the deduplicated storage from the original archive without overwriting differing inputs.

## Remaining inputs

The WGS read accession remains to be added. Commands and saved options document the runs; standalone job scripts are unavailable. Read alignment and mutation calling are not rerun by this repository. Single-cell inputs and count-matrix reconstruction are documented in ../single-cell/README.md.

## Software-version assumption

Version documentation assumes the latest public breseq release available on each recorded run date. The dates come from the original log entries, rather than file-copy timestamps. The [official release history](https://github.com/barricklab/breseq/releases) gives the following assignments:

| Recorded run dates | Runs | Assumed version | Release date |
|---|---:|---|---|
| 6–7 November 2023 and 2 January 2024 | 76 | [0.38.1](https://github.com/barricklab/breseq/releases/tag/v0.38.1) | 14 April 2023 |
| 16/24 April and 11 June 2024 | 30 | [0.38.3](https://github.com/barricklab/breseq/releases/tag/v0.38.3) | 4 February 2024 |
| 5 August and 22 October 2024 | 43 | [0.39.0](https://github.com/barricklab/breseq/releases/tag/v0.39.0) | 10 July 2024 |

Version 0.38.2 was released on 7 January 2024, after the January study runs. These date-based assignments are assumptions, not versions recovered from the executable.
