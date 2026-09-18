# Numerical results

36 primary IC50 comparisons; 369/369 acceptance checks passed. Source Data retain full precision.

IC50 comparisons use log10 measurements, culture-ID-aligned paired t-tests for Figures 2/3 and Supplementary Figures 3/6, and Welch tests where matching is unconfirmed. Holm correction is applied within panel and drug (three comparisons in Figure 2, six in Figure 3, one in each other tested panel). Ratios compare Group2 with Group1; confidence intervals are pointwise. MIC ratio intervals estimate uncertainty using the stated log-scale t method, without multiplicity adjustment; MIC has no hypothesis tests or significance calls. Figure 5a-d and Supplementary Figure 10a-c are technical replicates of one culture per strain.

The fit-perturbation sensitivity reports empirical 2.5th-97.5th percentiles from 4,000 culture/fit resamples (seed 20260915); primary inference uses the t-tests and t-based intervals. Other sensitivity p-values have explicitly identified adjustment families.

## Primary comparisons

| Panel | Drug | Comparison | Ratio (pointwise 95% CI) | Holm p |
|---|---|---|---|---|
| Figure 2b-d | Levofloxacin | PA/P | 0.782 (0.671-0.91) | p = 0.010 |
| Figure 2b-d | Levofloxacin | PC/P | 1.13 (0.837-1.53) | p = 0.375 |
| Figure 2b-d | Levofloxacin | PL/P | 36.1 (26.1-50.1) | p < 0.0001 |
| Figure 2b-d | Amikacin | PA/P | 12.2 (8.03-18.5) | p < 0.0001 |
| Figure 2b-d | Amikacin | PC/P | 1.1 (0.911-1.33) | p = 0.554 |
| Figure 2b-d | Amikacin | PL/P | 1 (0.744-1.34) | p = 0.999 |
| Figure 2b-d | Cefepime | PA/P | 0.79 (0.64-0.975) | p = 0.032 |
| Figure 2b-d | Cefepime | PC/P | 1.61 (1.13-2.28) | p = 0.026 |
| Figure 2b-d | Cefepime | PL/P | 2.08 (1.59-2.72) | p = 0.0005 |
| Figure 3b-d | Levofloxacin | PL/P | 36.1 (26.1-50.1) | p < 0.0001 |
| Figure 3b-d | Levofloxacin | PLA/P | 15.2 (10.6-21.9) | p < 0.0001 |
| Figure 3b-d | Levofloxacin | PLAC/P | 17.3 (10.7-28.2) | p < 0.0001 |
| Figure 3b-d | Levofloxacin | PLA/PL | 0.422 (0.294-0.605) | p = 0.001 |
| Figure 3b-d | Levofloxacin | PLAC/PL | 0.48 (0.297-0.776) | p = 0.014 |
| Figure 3b-d | Levofloxacin | PLAC/PLA | 1.14 (0.824-1.57) | p = 0.389 |
| Figure 3b-d | Amikacin | PL/P | 1 (0.744-1.34) | p = 0.999 |
| Figure 3b-d | Amikacin | PLA/P | 8.43 (6.63-10.7) | p < 0.0001 |
| Figure 3b-d | Amikacin | PLAC/P | 2.61 (1.58-4.31) | p = 0.006 |
| Figure 3b-d | Amikacin | PLA/PL | 8.43 (6.47-11) | p < 0.0001 |
| Figure 3b-d | Amikacin | PLAC/PL | 2.61 (1.49-4.58) | p = 0.008 |
| Figure 3b-d | Amikacin | PLAC/PLA | 0.31 (0.194-0.494) | p = 0.001 |
| Figure 3b-d | Cefepime | PL/P | 2.08 (1.59-2.72) | p = 0.001 |
| Figure 3b-d | Cefepime | PLA/P | 1.24 (0.902-1.71) | p = 0.319 |
| Figure 3b-d | Cefepime | PLAC/P | 1.39 (1-1.92) | p = 0.145 |
| Figure 3b-d | Cefepime | PLA/PL | 0.598 (0.431-0.83) | p = 0.025 |
| Figure 3b-d | Cefepime | PLAC/PL | 0.667 (0.537-0.827) | p = 0.011 |
| Figure 3b-d | Cefepime | PLAC/PLA | 1.12 (0.767-1.62) | p = 0.527 |
| Supplementary Figure 3b-d | Levofloxacin | PCr/P | 4 (2.04-7.82) | p = 0.003 |
| Supplementary Figure 3b-d | Amikacin | PCr/P | 1.77 (0.996-3.16) | p = 0.051 |
| Supplementary Figure 3b-d | Cefepime | PCr/P | 109 (13-922) | p = 0.002 |
| Supplementary Figure 4b-d | Levofloxacin | ATEC-C/ATEC | 1.09 (0.988-1.2) | p = 0.076 |
| Supplementary Figure 4b-d | Amikacin | ATEC-C/ATEC | 0.767 (0.502-1.17) | p = 0.192 |
| Supplementary Figure 4b-d | Cefepime | ATEC-C/ATEC | 0.811 (0.597-1.1) | p = 0.159 |
| Supplementary Figure 6a-c | Levofloxacin | P-unt/P | 1.22 (1.13-1.32) | p = 0.004 |
| Supplementary Figure 6a-c | Amikacin | P-unt/P | 0.711 (0.566-0.892) | p = 0.017 |
| Supplementary Figure 6a-c | Cefepime | P-unt/P | 0.677 (0.517-0.885) | p = 0.019 |

## Correction-family sensitivity

Figure 2 cefepime PA/P: primary Holm p = 0.032, per-drug Bonferroni p = 0.096, figure-wide Holm p = 0.128. Both comparisons depend on the correction family.
Figure 2 cefepime PC/P: primary Holm p = 0.026, per-drug Bonferroni p = 0.040, figure-wide Holm p = 0.066. Both comparisons depend on the correction family.

The Figure 3 cefepime PL-versus-PLA comparison has a primary adjusted p-value of 0.025. All adjusted values, including nonsignificant comparisons, are available in the tables.

Supplementary Figure 4 cefepime changes from the historical significant call to a nonsignificant log-scale comparison. Untreated passage changes IC50 in all three drugs, as quantified by the reported ratios.

## Mutant IC50, descriptive

| Strain | Drug | Ratio to MG | Technical range / MG geometric mean |
|---|---|---|---|
| gatA | Levofloxacin | 1.01 | 0.965-1.04 |
| glvC | Levofloxacin | 1.43 | 1.1-1.64 |
| hipA | Levofloxacin | 1.51 | 1.42-1.65 |
| selB | Levofloxacin | 2.38 | 1.84-2.93 |
| rpoZ | Levofloxacin | 0.95 | 0.795-1.17 |
| ftsH | Levofloxacin | 0.849 | 0.749-0.905 |
| fimE | Levofloxacin | 1.55 | 1.04-2.09 |
| gatA | Amikacin | 0.891 | 0.812-1.02 |
| glvC | Amikacin | 1.07 | 0.924-1.26 |
| hipA | Amikacin | 0.799 | 0.785-0.821 |
| selB | Amikacin | 3.43 | 2.99-3.73 |
| rpoZ | Amikacin | 0.704 | 0.602-0.891 |
| ftsH | Amikacin | 0.472 | 0.345-0.57 |
| fimE | Amikacin | 0.877 | 0.649-1.04 |
| gatA | Cefepime | 0.97 | 0.836-1.17 |
| glvC | Cefepime | 1.08 | 1.01-1.15 |
| hipA | Cefepime | 0.851 | 0.798-0.943 |
| selB | Cefepime | 1.09 | 1.03-1.14 |
| rpoZ | Cefepime | 0.833 | 0.714-0.933 |
| ftsH | Cefepime | 1.23 | 1.12-1.35 |
| fimE | Cefepime | 0.857 | 0.838-0.882 |

Figure 5a reports six technical wells per strain, with doubling times in minutes.

## Supplementary Figure 11, conditional normalization

The calculation assumes equal plated volumes across time. If the actual schedule was 10 microlitres at baseline and 20 microlitres later, absolute fractions would be half these values. Relative strain comparisons are unchanged under a common volume schedule. The workbook and written account specify inconsistent volume normalization.

| Strain | Time (h) | Median fraction | Ratio to hipA |
|---|---|---|---|
| MG hipA | 3 | 0.000131 | 1 |
| MG hipA | 7 | 1.09e-05 | 1 |
| MG hipA+ftsH | 3 | 0.00043 | 3.3 |
| MG hipA+ftsH | 7 | 3.38e-05 | 3.1 |
| MG hipA+rpoZ | 3 | 0.000805 | 6.16 |
| MG hipA+rpoZ | 7 | 2.58e-05 | 2.37 |
| MG hipA+selB | 3 | 0.000156 | 1.19 |
| MG hipA+selB | 7 | 6.31e-06 | 0.579 |
| MG | 3 | 1.59e-05 | 0.122 |
| MG | 7 | 8.8e-08 | 0.00808 |

The recorded zero is retained as censored (n = 4), with a separate plotted detection-limit marker. The WT 7 h median is 8.7971e-08; normal-scaled MAD ranges from 7.27104e-08 to 8.5948e-08 over the censored interval. The displayed MAD uses the one-colony upper bound. Missing Figure 3 counts are recorded separately from this observed zero.

Unresolved records and figure coverage are listed in the repository README and `docs/FIGURE_COVERAGE.md`.
