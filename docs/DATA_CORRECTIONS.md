# Data processing and source records

## Experimental units and survival

IC50 comparisons in Supplementary Figure 3 pair P1–P6 with PCr1–PCr6; Supplementary Figure 6 pairs P1–P4 with P-unt1–P-unt4. Each evolved or passaged culture is matched to its measured ancestor by culture number. P7–P10 remain in the Supplementary Figure 3 source data but have no evolved partners and are excluded from its comparison and plotted parent summary. Tests use log10 values and Holm adjustment within panel and antibiotic. Supplementary Figure 4 uses Welch tests because culture matching is not established.

For cefepime culture 10 at day 20, SurvivalData_26.xlsx (PC row 401) and bigcfu070523_26.xlsx (singleDrug row 401) record 16 colonies at dilution 10^-5 and 20 µL plated volume, giving 80,000,000 CFU/mL. The matched pre-treatment value is 150,000,000 CFU/mL, giving 53.333333% survival. The ten-culture mean is 34.7224000558%, with sample SD 30.5038182291%. The [source record](../data/source_corrections_2026-09-18.json) identifies the workbook hashes and values.

The final sequential-evolution and sequencing time point is day 82. Every strain–antibiotic combination in the six fitted datasets contains zero-dose growth controls in addition to media-only background wells.

## MDK counts, volumes and exclusions

Supplementary Figure 11 uses a single 10 µL plating at every measured time point. Equal volume cancels in normalized fractions and detection limits. Earlier MDK assays used two 10 µL platings and later assays one 20 µL plating, with recorded 40 µL concentrated-sample exceptions. Calculations use each record’s volume and dilution or concentration factor, including tenfold concentration where applicable.

Figure 3 source cells H62 and H146 represent zero-colony plates. These cells are blank in the original workbook; derived fields specify count 0, dilution exponent −1, factor 10 and 20 µL plating. One-colony fraction thresholds are 2.38095238e-8 and 2.22222222e-8. Zero and missing observations have separate status/count fields. Both zero records contribute to the primary summaries, whose parent medians and MADs are invariant over the censored intervals. Possible pellet loss for the second-day 3-hour record is unconfirmed. Omitting that record changes the median from 4.80769231e-7 to 5.00000000e-7 (4%) and the scaled MAD from 1.18402261e-7 to 8.23667899e-8.

ATEC1/2/4 experiment-2 measurements and the failed ATEC-C3 technical series were excluded and repeated in experiment 3. The experimental account attributes the ATEC-C3 failure to inoculation error; the available plate screenshot alone does not establish the cause. The exclusion lists and drug-specific filters are recorded in `pipeline.py`.

## Genomic tables and variant identity

Inputs comprise 149 mutation-call TSVs and 17 processed CSVs. The complete PbEc table contains 83 rows. The 17 processed exports agree with the notebook calculations before the annotation and allele transformations applied in later data-processing cells. Reference sequences, commands and saved options cover all 149 runs; date-based breseq version assumptions are documented in `data/genomics/README.md`.

Culture 2 endpoint frequencies are GyrA S83L, 1.000000; HipA, 0.154478; selB, 0.183824; fimB/fimE, 0.178571; ftsH, 0.187417; and gatA, 0.165929. The latter alleles are present at low frequency.

Plots identify variants by coordinate and allele. Distinct mutations remain separate even when they share a short label:

| Plot | Variant identities and source frequencies |
|---|---|
| S4f, culture 7, atpD | S342R at contig 39:25113, 0.43; indel at contig 39:25321, 0.27 |
| S4f, GlnW | Contig 22 positions 86211/86212: culture 6, 0.13/0.14; culture 7, 0.14/0.15 |
| S6d, culture 2, yhaC/rnpB | Positions 3269943/3270068: day 16, 0/0.10; day 82, 0.12/0 |
| S9, PLAC_01, cpxA | Positions 4104926/4104944, 0.40/0 |
| S9, PLAC_06, ftsH | Positions 3326144/3326602, 0.47/0; the insertion has its own label |
| S9, PLA_07, ampH/sbmA | Positions 396591/396624, 0.66/0.13 |

The [crosswalk](../data/genomics/plot_tables/panel_source_crosswalk.csv) records all 1,099 plotted population/variant entries, source labels, display aliases and coordinates. The GlnW sites remain separate because their raw frequencies differ and phasing is not established. Figure 4a labels the GyrA G81D allele; Figure 4c colors use the Figure 3 MDK phenotype classifications. Figure 4 and Supplementary Figure 8 trajectories use equally spaced sampled days, as stated in their captions.

## Single-cell analysis and enrichment

Count-matrix analysis retains 48,883 cells and nine clusters, producing 6,543 cluster DGE rows, 2,513 sample-versus-rest rows and 1,222/863 culture-versus-parent rows. Source cluster IDs 0–8 map to published clusters 1–9. Six samples are technical libraries, two from each of three culture samplings, separately probed and processed through microfluidics.

Supplementary Table 3 contains sample/cluster results. Table 4 contains culture comparisons and separately identified targeted marker checks. Supplementary Figures 12/13 use Excel chart annotations matched uniquely to full DGE rows, including display labels and colors. Figure 6 uses an absolute-log2-fold-change plotting cutoff of 0.5; its 91 shared-gene highlights use 0.6. Values below 10^-304 share a display floor.

rplJ has higher mean expression in the parent than in either evolved culture. In clusters 1 and 2, nominal P-values are 0.00626 and 0.00526, but both Bonferroni-adjusted P-values are 1. The hipA-positive fraction is higher in clusters 4 and 9 than in the remaining clusters. Mean expression is higher in cluster 9 and slightly lower in cluster 4 (log2 fold changes 0.597 and −0.136).

The default culture-7 DGE table excludes hipA because its detection frequency is below 1%. A targeted comparison without that filter gives log2 fold change 0.721 and adjusted P = 0.00377. Targeted checks retain Bonferroni correction over 21,701 assay features and are reported separately from the default DGE and plotted gene set.

EcoCyc text exports agree with all 746 overlapping workbook rows in term, P-value and matched genes. Supplementary Table 4 and Source Data include all 133 parent-up export rows. EcoCyc 29.0 and default annotation-specific backgrounds are assumptions based on the analysis dates and documented standard workflow; see the [enrichment record](../data/single-cell/enrichment/README.md).
