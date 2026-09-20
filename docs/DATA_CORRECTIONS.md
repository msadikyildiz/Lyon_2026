# Source corrections

## 18 September 2026

Adam Lyon confirmed that the parent IC50 measurements in Supplementary Figures 3 and 6 came from each lineage's unique ancestor. Comparisons now use paired log10 t-tests, aligned by culture number. Supplementary Figure 3 uses P1–P6 with PCr1–PCr6. P7–P10 remain in the raw cache and Source Data but have no corresponding evolved cultures in this panel, so are excluded from its comparison and plotted parent summary. Supplementary Figure 6 uses P1–P4 with P-unt1–P-unt4. Holm adjustment remains within panel and antibiotic. Supplementary Figure 4 continues to use Welch tests.

The source workbooks SurvivalData_26.xlsx (PC row 401) and bigcfu070523_26.xlsx (singleDrug row 401) restore the omitted day-20 post-treatment observation for cefepime culture 10: 16 colonies at dilution 10^-5, plated volume 20 µL, and 80,000,000 CFU/mL. The corresponding pre-treatment value is 150,000,000 CFU/mL, giving 53.333333% survival. The complete ten-culture day-20 mean is 34.7224000558%, with sample SD 30.5038182291%. Only the recovered row was added to the two canonical survival workbooks; their other sheets were preserved. The [source record](../data/source_corrections_2026-09-18.json) gives the supplied workbook hashes and recovered values.

The final sequential-evolution and sequencing time point is day 82. Fourteen terminal sample entries in the genomic metadata and the Figure 4/Supplementary Figure 8 notebook day mappings were corrected from 81 to 82. The retained Supplementary Figure 6 heatmap has the same terminal-label correction; its data pixels are unchanged. Amino-acid positions, including GyrA G81D, are unchanged.

Every strain–antibiotic combination in the six fitted datasets contains zero-dose wells. These are the no-antibiotic growth controls, in addition to media-only background wells. No dose-response measurements or fits changed.

The double-mutant median ratios and qualified SelB interpretation were confirmed. Supplementary Figure 11 plated volumes and concentration factors remain unresolved. Complete genomic exports and the Rosenthal group's single-cell data and code are still required, as listed in the README.

## 20 September 2026

The complete genomic input set is now included: 149 mutation-call TSVs and 17 processed CSVs supplied by Erdal Toprak. All 83 PbEc rows are regenerated from the calls, replacing the partial 50-row saved-display export in Source Data. The regenerated PLAC and combined-lineage tables also replace their previous dependency placeholders. Existing endpoint values agree with the regenerated tables, and all 17 processed exports agree before the notebooks' later annotation and allele corrections.

The Culture 2 endpoint calls confirm GyrA S83L at 1.000000, HipA at 0.154478, selB at 0.183824, fimB/fimE at 0.178571, ftsH at 0.187417 and gatA at 0.165929. These support describing the latter alleles as present at low frequency. Full reference-genome/caller provenance and sequencing accessions remain needed for reproduction from reads. Single-cell inputs and their subsequent reconstruction are documented in data/single-cell/README.md.

## Confirmed records and genomic plotting, 20 September 2026

The Supplementary Figure 11 plates all used a single 10 µL volume. This confirmation removes the provisional volume assumption without changing the normalized fractions. Recorded 40 µL concentrated-sample exceptions were intentional. Earlier MDK assays used two 10 µL platings; later assays used one 20 µL plating. Calculations retain each recorded schedule and the tenfold concentration where applicable.

Figure 3 cells H62 and H146 were confirmed zero-colony plates. The original workbook remains unchanged; derived fields record zero counts, dilution exponent −1, factor 10, and 20 µL plating. Their one-colony fraction thresholds are 2.38095238e-8 and 2.22222222e-8. Zero and missing observations have separate status/count fields. Both records remain in the primary summaries. Parent medians and MADs are invariant over the censored intervals. Omitting the possibly mishandled second-day 3-hour record raises the median from 4.80769231e-7 to 5.00000000e-7 (4%); the scaled MAD changes from 1.18402261e-7 to 8.23667899e-8. Pellet loss remains unconfirmed.

The ATEC-C3 experiment-2 triplicate was reported to have no bacteria and was repeated in experiment 3. ATEC1/2/4 experiment-2 errors were also repeated in experiment 3. This is the experimenter's clarification; the available screenshot alone does not establish the cause.

Genomic plotting now keys variants by coordinate and allele, preserving distinct sites under shared labels. The original mean-based pivot merged different mutations. Corrections include:

| Plot | Historical aggregation | Separate source frequencies |
|---|---|---|
| S4f, culture 7 | A gene-product alias labeled both atpD variants S342R and averaged them to 0.35 | S342R at contig 39:25113, 0.43; indel at contig 39:25321, 0.27 |
| S4f, GlnW | Adjacent positions 86211/86212 shared a row; culture-6/7 means were 0.135/0.145 | Culture 6: 0.13/0.14; culture 7: 0.14/0.15; contig 22 |
| S6d, culture 2 | Two yhaC/rnpB positions averaged to 0.05 at day 16 and 0.06 at day 82 | Day 16: 0/0.10; day 82: 0.12/0; positions 3269943/3270068 |
| S9, PLAC_01 | cpxA variants averaged to 0.20 | 0.40 and 0, positions 4104926/4104944 |
| S9, PLAC_06 | Distinct ftsH variants averaged to 0.235 | 0.47 and 0, positions 3326144/3326602; the insertion retains its own label |
| S9, PLA_07 | ampH/sbmA variants averaged to 0.395 | 0.66 and 0.13, positions 396591/396624 |

The full [crosswalk](../data/genomics/plot_tables/panel_source_crosswalk.csv) records all 1,099 plotted population/variant entries, original labels, displayed aliases and coordinates. Different sites found only in different cultures also receive distinct rows. The GlnW sites remain separate because raw frequencies differ and phasing is not established by these tables. Source table inputs and historical matrix calculations remain available. Figure 4a now labels GyrA G81D correctly; Figure 4c phenotype colors follow the existing Figure 3 MDK classifications. Figure 4d retains the supplied schematic pending the final export. Trajectory plots retain the original notebook convention of equally spaced sampled days; Figure 4 and Supplementary Figure 8 captions state this convention.

## Single-cell reconstruction and interpretation, 20 September 2026

The supplied matrices and original Rmd reproduce all 48,883 retained cells, all nine cluster sizes, 6,543 cluster DGE rows, 2,513 sample-versus-rest rows, and 1,222/863 culture-versus-parent rows. The source cluster IDs 0–8 are shifted to the published 1–9. Six samples are technical libraries, two from each of three culture samplings, separately probed and processed through microfluidics.

The full supplied tables replace missing-data entries in Source Data. Supplementary Table 3 now contains sample/cluster results, while Table 4 contains culture comparisons and separately identified targeted marker checks. The original workbook remains unchanged. S12/13 point annotations are read directly from the Excel charts and matched uniquely to the recomputed DGE rows, preserving manually assigned display labels and colors. The default significance/absolute-log2-fold-change thresholds are unchanged. Figure 6 legend and caption distinguish the 0.5 plotting cutoff from the 0.6 cutoff that exactly reproduces the original 91 shared-gene highlights. Values below 10^-304 share a display floor.

The reconstruction agrees with the supplied expression tables and clarifies the manuscript's statistical wording. rplJ has higher mean expression in the parent than in either evolved culture. In clusters 1 and 2, its nominal P-values are 0.00626 and 0.00526, but both Bonferroni-adjusted P-values are 1. The hipA-positive cell fraction is higher in both clusters 4 and 9 than in the remaining clusters. Mean expression is higher in cluster 9 and slightly lower in cluster 4 (log2 fold changes 0.597 and −0.136). The original qualitative distribution statement is therefore consistent with detection frequency; a negative mean fold change alone does not contradict it.

The culture-7 hipA omission was caused by the 1% detection-frequency filter. A targeted check without that filter gives log2 fold change 0.721 and adjusted P = 0.00377. These checks retain the supplied correction over 21,701 assay features and are reported separately from the default DGE and plotted gene set. An independent calculation directly from H5 counts reproduces the selected fold changes and P-values. No source counts, clustering or default DGE calculations were changed.

Enrichment worksheets include full result tables and gene memberships. Their generating tool/settings, database version and tested background are still required. The public GSE314756 endpoint reported a 21 December 2026 release date on 20 September 2026, despite the accompanying description that the counts were public. No reviewer credential is published.
