# Source corrections

## 18 September 2026

Adam Lyon confirmed that the parent IC50 measurements in Supplementary Figures 3 and 6 came from each lineage's unique ancestor. Comparisons now use paired log10 t-tests, aligned by culture number. Supplementary Figure 3 uses P1–P6 with PCr1–PCr6. P7–P10 remain in the raw cache and Source Data but have no corresponding evolved cultures in this panel, so are excluded from its comparison and plotted parent summary. Supplementary Figure 6 uses P1–P4 with P-unt1–P-unt4. Holm adjustment remains within panel and antibiotic. Supplementary Figure 4 continues to use Welch tests.

The source workbooks SurvivalData_26.xlsx (PC row 401) and bigcfu070523_26.xlsx (singleDrug row 401) restore the omitted day-20 post-treatment observation for cefepime culture 10: 16 colonies at dilution 10^-5, plated volume 20 µL, and 80,000,000 CFU/mL. The corresponding pre-treatment value is 150,000,000 CFU/mL, giving 53.333333% survival. The complete ten-culture day-20 mean is 34.7224000558%, with sample SD 30.5038182291%. Only the recovered row was added to the two canonical survival workbooks; their other sheets were preserved. The [source record](../data/source_corrections_2026-09-18.json) gives the supplied workbook hashes and recovered values.

The final sequential-evolution and sequencing time point is day 82. Fourteen terminal sample entries in the genomic metadata and the Figure 4/Supplementary Figure 8 notebook day mappings were corrected from 81 to 82. The retained Supplementary Figure 6 heatmap has the same terminal-label correction; its data pixels are unchanged. Amino-acid positions, including GyrA G81D, are unchanged.

Every strain–antibiotic combination in the six fitted datasets contains zero-dose wells. These are the no-antibiotic growth controls, in addition to media-only background wells. No dose-response measurements or fits changed.

The double-mutant median ratios and qualified SelB interpretation were confirmed. Supplementary Figure 11 plated volumes and concentration factors remain unresolved. Complete genomic exports and the Rosenthal group's single-cell data and code are still required, as listed in the README.

## 20 September 2026

The complete genomic input set is now included: 149 mutation-call TSVs and 17 processed CSVs supplied by Erdal Toprak. All 83 PbEc rows are regenerated from the calls, replacing the partial 50-row saved-display export in Source Data. The regenerated PLAC and combined-lineage tables also replace their previous dependency placeholders. Existing endpoint values agree with the regenerated tables, and all 17 processed exports agree before the notebooks' later annotation and allele corrections.

The Culture 2 endpoint calls confirm GyrA S83L at 1.000000, HipA at 0.154478, selB at 0.183824, fimB/fimE at 0.178571, ftsH at 0.187417 and gatA at 0.165929. These support describing the latter alleles as present at low frequency. Full reference-genome/caller provenance and sequencing accessions remain needed for reproduction from reads. Single-cell inputs remain separate and outstanding.
