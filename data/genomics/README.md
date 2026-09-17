# Genomic inputs

`reference/metadata_complete.csv` is the available sample metadata. `recovered_tables/index.json` maps each saved notebook table to a portable JSON table with row counts, completeness and notebook cell provenance. Eight entries contain complete displayed tables, including repeated copies for the extended heatmaps; one contains only 50 of 83 rows. Displayed endpoint frequencies retain the precision printed by the original notebooks. Long annotation strings or frequency vectors may have been truncated in that display, so these exports cannot replace the original sequencing tables.

The original plotting and processing code is retained in the notebooks listed in `notebook_inventory.json`. The first cell locates the repository root so paths work when a notebook is opened in its figure directory. The automated reproduction command uses recovered tables for Source Data and retains the existing genomic image panels in figure assemblies.

## Input requirements

Genomic figure regeneration requires either the per-sample mutation tables below or the full processed frequency tables used by the plotting notebooks, with sample/culture/day and mutation identities. The processed tables should include full-precision frequencies and the full traced-allele vectors, particularly PLAC for Figure 4 and Supplementary Figures 8/9, and the missing 33 rows for the PbEc table in Supplementary Figure 4.

The original notebook input layout is:

```text
data/genomics/reference/metadata_complete.csv
data/genomics/out/<FolderDate>/<source_file>/output/output.gd.tsv
```

`FolderDate` and `source_file` are metadata columns. The tables require mutation-calling/reference provenance and the corresponding plotting notebook version. If rerunning the mutation calls is intended, reference genome, caller version/options and read accessions are also needed. The supplied notebooks cover plotting and downstream table processing, not a complete raw-read alignment/calling workflow.

Single-cell data are separate: Figure 6 and Supplementary Figures 12/13 need their expression/metadata objects, embeddings, differential-expression/enrichment tables and generating code from the Rosenthal analysis. None of those inputs is present in this repository.
