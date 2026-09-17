# Input provenance

The original dose-response and plate-layout workbooks are in `working/analysis/stats-rework/biohpc-pull/data`. Their original trailing OD column is preserved because the historical preprocessing removes that column explicitly. The six fitted caches retain the culture IDs, fit grids, bootstrap draws and dose-level records used in the reviewed analysis. plategig source is vendored at commit `88839a2`; its MIT license is included.

Growth, survival and MDK workbooks are retained under their figure directories. The three supplementary manuscript tables were exported to `data/supplementary_tables.json`; no manuscript DOCX is needed to regenerate Source Data. Existing figure panels were extracted without altering their pixels to `data/figure-assets`, with their original hashes and panel identities recorded in `assembly_assets.json`.

Genomic metadata and saved displayed tables are described in `data/genomics/README.md`. Notebook outputs and execution metadata were cleared; scientific code is retained with portable paths and a repository-root bootstrap. Five survival notebooks select the bundled DejaVu Sans font explicitly, matching the original Nimbus Roman fallback on the validation machine. The original per-sample mutation-call files and complete single-cell analysis inputs remain unavailable.

`data/INPUT_SHA256SUMS` records the supplied data, reference notebooks, cache files and vendored source. Verify from the repository root with `shasum -a 256 -c data/INPUT_SHA256SUMS`. Generated results and runtime logs are excluded from that input checksum list.
