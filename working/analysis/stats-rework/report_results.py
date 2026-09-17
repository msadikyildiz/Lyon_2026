"""Generate the numerical review note directly from full-precision outputs."""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent


def p_text(value):
    if value < .0001:
        return 'p < 0.0001'
    return f'p = {value:.4f}' if value < .001 else f'p = {value:.3f}'


def main():
    t = pd.read_csv(HERE / 'out/final_statistics.csv')
    ic = t[t.value == 'IC50']
    f5 = pd.read_csv(HERE / 'out/fig5_descriptive.csv').query('value == "IC50"')
    m = pd.read_csv(HERE / 'out/supp11_summary.csv')
    checks = pd.read_csv(HERE / 'out/acceptance_checks.csv')
    lines = ['# Numerical results', '',
             f'{len(ic)} primary IC50 comparisons; {len(checks)}/{len(checks)} acceptance checks passed. '
             'Source Data retain full precision.', '',
             'IC50 comparisons use log10 measurements, culture-ID-aligned paired t-tests for Figures 2 and 3, '
             'and Welch tests where matching is unconfirmed. Holm correction is applied within panel and drug '
             '(three comparisons in Figure 2, six in Figure 3, one in each other tested panel). '
             'Ratios compare Group2 with Group1; confidence intervals are pointwise. MIC ratio intervals '
             'estimate uncertainty using the stated log-scale t method, without multiplicity adjustment; '
             'MIC has no hypothesis tests or significance calls. '
             'Figure 5a-d and Supplementary Figure 10a-c are technical replicates of one culture per strain.', '',
             'The fit-perturbation sensitivity reports empirical 2.5th-97.5th percentiles from 4,000 '
             'culture/fit resamples (seed 20260915); primary inference uses the t-tests and t-based intervals. '
             'Other sensitivity p-values have explicitly identified adjustment families.', '',
             '## Primary comparisons', '',
             '| Panel | Drug | Comparison | Ratio (pointwise 95% CI) | Holm p |',
             '|---|---|---|---|---|']
    for _, r in ic.iterrows():
        lines.append(f'| {r.Panel} | {r.Antibiotic} | {r.Group2}/{r.Group1} | '
                     f'{r.ratio:.3g} ({r.ratio_lo:.3g}-{r.ratio_hi:.3g}) | {p_text(r.p_holm)} |')
    lines += ['', '## Correction-family sensitivity', '']
    for g in ['PA', 'PC']:
        r = ic[(ic.Panel == 'Figure 2b-d') & (ic.Antibiotic == 'Cefepime') & (ic.Group2 == g)].iloc[0]
        lines.append(f'Figure 2 cefepime {g}/P: primary Holm {p_text(r.p_holm)}, '
                     f'per-drug Bonferroni {p_text(r.p_bonferroni)}, '
                     f'figure-wide Holm {p_text(r.p_holm_figure)}. Both comparisons depend on the correction family.')
    lines += ['', 'Reviewer 2 identifies Figure 3 cefepime PL versus PLA. Its primary adjusted '
              'p-value is 0.025, not the PL-versus-PLAC value. All adjusted values, including '
              'nonsignificant comparisons, are available in the tables.', '',
              'Supplementary Figure 4 cefepime changes from the historical significant call to '
              'a nonsignificant log-scale comparison. Untreated passage changes IC50 in all three '
              'drugs; describe the ratios rather than calling that control unchanged.', '',
              '## Mutant IC50, descriptive', '',
              '| Strain | Drug | Ratio to MG | Technical range / MG geometric mean |', '|---|---|---|---|']
    for _, r in f5[f5.Strain != 'MG'].iterrows():
        lines.append(f'| {r.Strain} | {r.Antibiotic} | {r.ratio_vs_MG:.3g} | {r.ratio_min:.3g}-{r.ratio_max:.3g} |')
    lines += ['', 'Figure 5a reports six technical wells per strain, with doubling times in minutes.', '',
              '## Supplementary Figure 11, conditional normalization', '',
              'The calculation assumes equal plated volumes across time. If the actual schedule '
              'was 10 microlitres at baseline and 20 microlitres later, absolute fractions would '
              'be half these values. Relative strain comparisons are unchanged under a common '
              'volume schedule. Adam must reconcile the workbook and written account.', '',
              '| Strain | Time (h) | Median fraction | Ratio to hipA |', '|---|---|---|---|']
    for _, r in m[m.Time > 0].iterrows():
        lines.append(f'| {r.strain} | {r.Time} | {r["median"]:.3g} | {r.ratio_to_hipA_median:.3g} |')
    r = m[m.n_below_detection > 0].iloc[0]
    lines += ['', f'The recorded zero is retained as censored (n = 4), with a separate plotted '
              f'detection-limit marker. The WT 7 h median is {r["median"]:.6g}; normal-scaled MAD '
              f'ranges from {r.mad_scaled_min:.6g} to {r.mad_scaled_max:.6g} over the censored interval. '
              'The displayed MAD uses the one-colony upper bound. Missing Figure 3 counts are '
              'recorded separately from this observed zero.', '',
              'Unresolved records and figure coverage are listed in the repository README and '
              '`docs/FIGURE_COVERAGE.md`.', '']
    (HERE / 'RESULTS_v2.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
