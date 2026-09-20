# Genomic inputs and calculations

The repository includes all 149 per-sample mutation-call TSVs referenced by `reference/metadata_complete.csv`, plus seven ATEC and ten PLAC processed CSV exports. These files were supplied by Erdal Toprak on 20 September 2026 in `Genomic_input_files.zip`. `input_manifest.json` records the archive checksum and each file's size, row count and SHA-256 checksum. Supplied files are preserved unchanged.

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

Each table is exported as CSV and typed JSON. `generated_tables/index.json` maps all 12 figure notebooks to their output tables and records executed cell numbers and notebook hashes. The 17 supplied processed exports match the rerun calculations before the notebooks' later annotation/allele corrections. Every available saved endpoint row also matches after accounting for truncated display text, empty-field formatting and numerical serialization. `generated_tables/validation.json` records these checks. The full 83-row PbEc table replaces the earlier 50-row display recovery in Source Data.

The original input layout is retained:

```text
data/genomics/reference/metadata_complete.csv
data/genomics/out/<FolderDate>/<source_file>/output/output.gd.tsv
data/genomics/data/processed/traced_alleles/<lineage>/<population>.csv
```

`recovered_tables/` preserves the historical notebook-display exports for comparison. They are no longer the Source Data input. The genomic plotting notebooks are included, but the automated figure assemblies still use their established image panels. Connecting regenerated genomic plots to those assemblies remains a separate task in the repository TODO table.

## Remaining inputs

Reproducing the upstream mutation calls requires reference-genome files/version, caller version/options and sequencing-read accessions. These are not supplied by this archive. Single-cell expression objects, metadata, embeddings, complete differential-expression/enrichment tables and generating code are separate inputs still needed for Figure 6 and Supplementary Figures 12/13.
