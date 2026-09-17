"""Export MDK records without converting missing colony counts to observed zeros."""
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.stats import median_abs_deviation
from manifest import HERE, ROOT, MDK


def read_records(figure):
    panel, rel, stem, note = MDK[figure]
    path = ROOT / 'working/figures' / rel
    d = pd.read_excel(path, sheet_name='Sheet1')
    d = d.rename(columns={'COUNT': 'Count', 'DILUTION': 'Dilution'})
    d['source_excel_row'] = d.index + 2
    d = d[d.Name.notna()].copy()
    d['count_status'] = np.where(d.Count.isna(), 'missing; requires confirmation',
                                  np.where(d.Count == 0, 'recorded zero; below detection', 'measured'))
    d['fraction_workbook'] = d.Fraction
    d['fraction_observed'] = d.Fraction.where(d.Count.notna() & (d.Count > 0))
    d['included_in_original_panel'] = True
    if figure == 'Figure 2':
        d.loc[d.Day == 1, 'included_in_original_panel'] = False
    d['replicate_original'] = d.Replicate
    if figure in ['Figure 2', 'Figure 3']:
        day = 3 if figure == 'Figure 2' else 2
        d.loc[(d.Day == day) & (d.Name == 'P1'), 'Replicate'] += 3
    ws = load_workbook(path, data_only=False, read_only=True)['Sheet1']
    d['cfu_formula'] = [ws[f'K{r}'].value for r in d.source_excel_row]
    d['formula_volume_mL'] = np.where((d.Count > 0) & (d.CFU > 0), d.Count / d.Factor / d.CFU, np.nan)
    d['record_note'] = note
    d['source_workbook'] = str(path.relative_to(ROOT))
    d['figure'] = figure
    return d


def main():
    for figure, (_, _, stem, _) in MDK.items():
        if figure == 'Supplementary Figure 11':
            continue
        d = read_records(figure)
        selected = d[d.included_in_original_panel]
        summaries = []
        for (name, time), g in selected.groupby(['Name', 'Time']):
            v = g.fraction_observed.dropna()
            summaries.append(dict(Name=name, Time=time, n_records=len(g), n_observed=len(v),
                                  n_missing=int(g.Count.isna().sum()), median=v.median(),
                                  mad_scaled=median_abs_deviation(v, scale='normal'),
                                  median_workbook=g.Fraction.median(),
                                  mad_scaled_workbook=median_abs_deviation(g.Fraction, scale='normal'),
                                  note='Observed-only summary is conditional where counts are missing'))
        s = pd.DataFrame(summaries)
        if figure in ['Figure 5', 'Supplementary Figure 10']:
            if (s.n_records != 3).any() or s.n_missing.any():
                raise ValueError(f'{figure}: expected three observed biological replicates')
            if not np.allclose(d.formula_volume_mL, .02, rtol=1e-12, atol=0):
                raise ValueError(f'{figure}: unexpected volume conversion')
        d.to_csv(HERE / 'out' / f'{stem}_records.csv', index=False)
        s.to_csv(HERE / 'out' / f'{stem}_summary.csv', index=False)
        print(f'{figure}: {len(d)} records; {selected.Count.isna().sum()} included missing counts')


if __name__ == '__main__':
    main()
