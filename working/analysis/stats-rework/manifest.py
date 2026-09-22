"""Panel identities and experimental units for the figure analyses.

Paths are relative to this repository. Panel letters follow the manuscript
legends, and are recorded in panel_manifest.csv.
"""
from pathlib import Path
import pandas as pd
from stats import split_strain

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'out'
DRUGS = ['Amikacin', 'Levofloxacin', 'Cefepime']
EXPECTED_IDS = {
    'fig2_paplpc': {g: list(map(str, range(1, 11))) for g in ['P', 'PA', 'PL', 'PC']},
    'fig3_plac': {g: list(map(str, range(1, 11))) for g in ['P', 'PL', 'PLA', 'PLAC']},
    'supp3_pcr': {'P': list(map(str, range(1, 11))), 'PCr': list(map(str, range(1, 7)))},
    'supp4_atec': {'ATEC': list(map(str, range(1, 7))),
                   'ATEC-C': list(map(str, range(1, 7))), 'ATEC-C-r': ['1']},
    'supp6_unt': {'P': list(map(str, range(1, 5))), 'P-unt': list(map(str, range(1, 5)))},
    'fig5_mutants': {g: ['1', '2', '3'] for g in
                     ['MG', 'gata', 'glvc', 'hipa', 'selb', 'rpoz', 'ftsh', 'fime']},
}


def comparison_cohort(dataset, table):
    """Use the six confirmed PCr parent pairs; preserve all rows in the raw cache."""
    if dataset != 'supp3_pcr':
        return table
    expected = {str(i) for i in range(1, 7)}
    for drug, rows in table.groupby('Antibiotic'):
        for group in ('P', 'PCr'):
            ids = rows.loc[rows['group'] == group, 'culturenumber'].astype(str)
            if ids.duplicated().any() or not expected.issubset(set(ids)):
                raise ValueError(f'{dataset}/{drug}/{group}: missing or duplicate confirmed pairs')
        if set(rows.loc[rows['group'] == 'PCr', 'culturenumber'].astype(str)) != expected:
            raise ValueError(f'{dataset}/{drug}: unexpected evolved lineage IDs')
    keep = (table['group'] != 'P') | table['culturenumber'].astype(str).isin(expected)
    return table.loc[keep].copy()

# Figure, letters, dataset, measure, render stem, drug order, statistics panel.
FIT_PANELS = [
    ('Figure 2', 'bcd', 'fig2_paplpc', 'IC50', 'Fig2', DRUGS, 'Figure 2b-d'),
    ('Figure 3', 'bcd', 'fig3_plac', 'IC50', 'Fig3',
     ['Levofloxacin', 'Amikacin', 'Cefepime'], 'Figure 3b-d'),
    ('Figure 5', 'bcd', 'fig5_mutants', 'IC50', 'Fig5', DRUGS, ''),
    ('Supplementary Figure 2', 'aaa', 'fig2_paplpc', 'IC50', 'Fig2', DRUGS, 'Figure 2b-d'),
    ('Supplementary Figure 2', 'bbb', 'fig2_paplpc', 'MIC', 'SuppFig2', DRUGS, 'Figure 2b-d'),
    ('Supplementary Figure 3', 'bcd', 'supp3_pcr', 'IC50', 'SuppFig3', DRUGS, 'Supplementary Figure 3b-d'),
    ('Supplementary Figure 4', 'bcd', 'supp4_atec', 'IC50', 'SuppFig4', DRUGS, 'Supplementary Figure 4b-d'),
    ('Supplementary Figure 6', 'abc', 'supp6_unt', 'IC50', 'SuppFig6',
     ['Amikacin', 'Cefepime', 'Levofloxacin'], 'Supplementary Figure 6a-c'),
    ('Supplementary Figure 7', 'aaa', 'fig3_plac', 'IC50', 'Fig3', DRUGS, 'Figure 3b-d'),
    ('Supplementary Figure 7', 'bbb', 'fig3_plac', 'MIC', 'SuppFig7', DRUGS, 'Figure 3b-d'),
    ('Supplementary Figure 10', 'abc', 'fig5_mutants', 'MIC', 'SuppFig10', DRUGS, ''),
]
MDK = {
    'Figure 2': ('e', 'Figure 2/E - MDK/mdk.P.PC.cultures-3days.xlsx', 'Fig2e_MDK',
                 'Record-specific volumes include 40 uL concentrated samples; analysis includes days 2 and 3'),
    'Figure 3': ('e', 'Figure 3/E - MDK/mdk.P.PLAC.cultures-2days.xlsx', 'Fig3e_MDK',
                 'Record-specific volumes; two zero-colony counts; possible handling loss at day 2, 3 h is unconfirmed'),
    'Figure 5': ('e', 'Figure 5 - mutants/E - MDK/MDK_mutants_correct.xlsx', 'Fig5e_MDK', ''),
    'Supplementary Figure 4': ('e', 'Supplemental Figure 4 - Pb/E - MDK/mdk.ATECc.cultures.xlsx',
                               'SuppFig4e_MDK', 'Recorded concentrated-sample plated volumes include 40 uL'),
    'Supplementary Figure 10': ('d', 'Supplemental Figure 10 - mutant prs and other/D - prs_hipA mdk/MDK_prs_mutant.xlsx',
                                'SuppFig10d_MDK', ''),
    'Supplementary Figure 11': ('all', 'Supplemental Figure 11 - double mutants/MDK_double_mutant.xlsx',
                                'SuppFig11_MDK_double_mutants', '10 uL at every measured time point; one censored count; 1/2 h not measured, 5 h excluded for dilution error'),
}


def panel_manifest():
    rows = []
    def add(fig, panel, measure, source, **kw):
        rows.append(dict(figure=fig, panel=panel, measure=measure, source=source, **kw))
    for fig, letters, ds, value, stem, drugs, stat_panel in FIT_PANELS:
        source = f'working/analysis/stats-rework/cache/{ds}.pkl'
        raw_suffix = '_mutants' if ds == 'fig5_mutants' else '_unt' if ds == 'supp6_unt' else ''
        for letter, drug in zip(letters, drugs):
            add(fig, letter, value, source, drug=drug, dataset=ds, units='ug/mL',
                output=f'out/figures-final/{stem}_{drug}_{value}', statistics_panel=stat_panel,
                raw_data=f'biohpc-pull/data/ODFinal{raw_suffix}.xlsx',
                unit=('technical dose-response series; biological n=1 per strain'
                      if ds == 'fig5_mutants' else 'culture'),
                status='standalone plot', note='Shared IC50 panel reuses primary adjustment family'
                if fig in ['Supplementary Figure 2', 'Supplementary Figure 7'] and value == 'IC50' else '')
    for fig, (panel, path, stem, note) in MDK.items():
        add(fig, panel, 'MDK survivor fraction', 'working/figures/' + path,
            units='fraction', output='out/figures-final/' + stem, drug='Cefepime',
            unit='biological culture; day and replicate retained',
            status='standalone plot' if fig in ['Figure 5', 'Supplementary Figure 10', 'Supplementary Figure 11'] else 'composite plot', note=note)
    add('Figure 5', 'a', 'doubling time',
        'working/figures/Figure 5 - mutants/A - doubling times/growth_022525_mutants.xlsx',
        units='min', output='out/figures-final/Fig5a_doubling_time',
        unit='six technical wells; biological n=1 per strain', status='standalone plot')
    for fig, panels, path, measure, units in [
        ('Figure 1', 'def', 'Figure 1/D-F/SurvivalData.xlsx', 'survival trajectory', 'percent'),
        ('Figure 3', 'a', 'Figure 3/A - survival/SurvivalData.xlsx', 'survival trajectory', 'percent'),
        ('Supplementary Figure 4', 'a', 'Supplemental Figure 4 - Pb/A - Survival/SurvivalData.xlsx', 'survival trajectory', 'percent'),
        ('Supplementary Figure 3', 'a', 'Supplemental Figure 3 - CefR/A - ResistanceEvo/pcr_evolution.xlsx', 'resistance trajectory', 'ug/mL'),
    ]:
        for panel in panels:
            add(fig, panel, measure, 'working/figures/' + path, units=units,
                unit='culture tracked over time', status='generated from data')
    notebooks = [
        ('Figure 2', 'f', 'Figure 2/F-H heatmaps/05-PA.ipynb'),
        ('Figure 2', 'g', 'Figure 2/F-H heatmaps/07-PL.ipynb'),
        ('Figure 2', 'h', 'Figure 2/F-H heatmaps/06-PC.ipynb'),
        ('Figure 4', 'abc', 'Figure 4/02-PLAC.ipynb'),
        ('Supplementary Figure 3', 'e', 'Supplemental Figure 3 - CefR/E - genetics/08-PCr.ipynb'),
        ('Supplementary Figure 4', 'f', 'Supplemental Figure 4 - Pb/F - Mutations/10-ATEC12.ipynb'),
        ('Supplementary Figure 5', 'a', 'Supplemental Figure 5 - Extended heatmaps/A/05-PA.ipynb'),
        ('Supplementary Figure 5', 'b', 'Supplemental Figure 5 - Extended heatmaps/B/07-PL.ipynb'),
        ('Supplementary Figure 5', 'c', 'Supplemental Figure 5 - Extended heatmaps/C/06-PC.ipynb'),
        ('Supplementary Figure 6', 'd', 'Supplemental Figure 6/B - heatmap/05-Punt.ipynb'),
        ('Supplementary Figure 8', 'abcde', 'Supplemental Figure 8 - extended traced alleles/02-PLAC.ipynb'),
        ('Supplementary Figure 9', 'all', 'Supplemental Figure 9 - big heatmap/13-bigheatmap.ipynb'),
    ]
    for fig, panels, path in notebooks:
        for panel in ([panels] if panels == 'all' else panels):
            add(fig, panel, 'mutation frequency', 'working/figures/' + path,
                units='fraction', unit='population sequencing sample',
                status='table and plot regenerated from mutation calls',
                note='Supp 5b: MG_LEV; Supp 5c: MG_CEF' if fig == 'Supplementary Figure 5' else '')
    for fig, panels, source in [
        ('Figure 1', 'abc', 'working/figures/Figure 1/A-C/Figure 1 A-C.pptx'),
        ('Figure 4', 'd', 'data/figure-assets/Figure4_schematic.json'),
        ('Supplementary Figure 1', 'all', 'working/figures/Supplemental Figure 1/cartoon - survival mechanisms.pptx'),
    ]:
        add(fig, panels, 'schematic', source, status='retain schematic', unit='not applicable', units='not applicable')
    add('Figure 2', 'a', 'representative dose-response curve',
        'working/figures/Figure 2/A-D - MIC/B-D/241011_Adam_mic_PAPLPC.ipynb',
        status='PA5 levofloxacin; notebook cell 29', units='OD and ug/mL', unit='dose-response series')
    for fig in ['Figure 6', 'Supplementary Figure 12', 'Supplementary Figure 13']:
        add(fig, 'all', 'single-cell RNA analysis', 'data/single-cell/input_manifest.json',
            status='count-matrix reconstruction and supplied-table validation', unit='cell; two technical libraries per culture')
    for i in [1, 2, 3, 4]:
        add(f'Supplementary Table {i}', 'all', 'table',
            'data/supplementary_tables.json' if i<3 else 'data/single-cell/original-tables/Single Cell Analysis Supp Tables.xlsx',
            status='original table in Source Data' if i<3 else 'complete expression/enrichment tables; EcoCyc assumptions in data/single-cell/enrichment/README.md')
    d = pd.DataFrame(rows).fillna('')
    if d.duplicated(['figure', 'panel', 'measure', 'drug']).any():
        raise ValueError('Duplicate panel identity')
    for p in d.source.unique():
        if not (ROOT / p).exists():
            raise FileNotFoundError(p)
    return d


def experimental_units():
    rows = []
    for ds in EXPECTED_IDS:
        gf = split_strain(pd.read_pickle(HERE / 'cache' / f'{ds}.pkl'))
        raw = pd.read_pickle(HERE / 'cache' / f'{ds}__df_analysis.pkl')
        for _, r in gf.iterrows():
            sub = raw[(raw.Strain == r.Strain) & (raw.Antibiotic == r.Antibiotic)]
            row = dict(dataset=ds, strain=r.Strain, group=r['group'], culture_id=str(r.culturenumber),
                       drug=r.Antibiotic, IC50=r.IC50, MIC=r.MIC,
                       biological_unit=(r['group'] if ds == 'fig5_mutants' else r.Strain),
                       unit_type='technical series' if ds == 'fig5_mutants' else 'culture',
                       lineage_id=str(r.culturenumber) if ds in ('fig2_paplpc', 'fig3_plac', 'supp3_pcr', 'supp6_unt') else '',
                       matching=('same-numbered ancestor and descendant cultures' if ds in ('supp3_pcr', 'supp6_unt')
                                 else 'same-numbered longitudinal cultures' if ds.startswith(('fig2_', 'fig3_')) else 'not established'),
                       included_in_comparison=not (ds == 'supp3_pcr' and r['group'] == 'P' and int(r.culturenumber) > 6),
                       analysis_source=f'cache/{ds}__df_analysis.pkl')
            for field in ['Experiment', 'Day', 'Plate', 'Plate_ID', 'Well', 'Row', 'Column']:
                row[field] = '|'.join(map(str, sub[field].dropna().unique())) if field in sub else ''
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    panels = panel_manifest()
    units = experimental_units()
    OUT.mkdir(exist_ok=True)
    panels.to_csv(OUT / 'panel_manifest.csv', index=False)
    units.to_csv(OUT / 'experimental_units.csv', index=False)
    print(f'{len(panels)} panel blocks; {len(units)} fitted experimental units')
