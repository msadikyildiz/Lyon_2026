"""Recalculate Figure 1f from paired pre/post counts and render its print panel."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', str(Path.home() / '.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'out'


def records():
    raw = pd.read_excel(ROOT / 'working/figures/Figure 1/D-F/SurvivalData.xlsx', sheet_name='PC')
    d = raw[(raw.Strain == 'MG1655') & raw.Culture.isin(range(1, 11)) & raw.Day.between(1, 20.5)].copy()
    d['culture'] = d.Culture.astype(int)
    d['day'] = np.floor(d.Day).astype(int)
    d['phase'] = np.where(d.Day == d['day'], 'before', 'after')
    if d.duplicated(['culture', 'day', 'phase']).any():
        raise ValueError('Duplicate survival observations')
    pairs = d.pivot(index=['culture', 'day'], columns='phase', values='CFU')
    expected = pd.MultiIndex.from_product([range(1, 11), range(1, 21)], names=['culture', 'day'])
    if set(pairs.index) != set(expected) or pairs.isna().any().any():
        raise ValueError('Missing pre/post culture-day pair')
    if (pairs['before'] <= 0).any() or (pairs['after'] < 0).any():
        raise ValueError('Invalid colony concentration')
    pairs = pairs.reindex(expected).reset_index().rename(columns={'before': 'pre_CFU_mL', 'after': 'post_CFU_mL'})
    pairs['survival_percent'] = 100 * pairs.post_CFU_mL / pairs.pre_CFU_mL
    return pairs


def main():
    d = records()
    summary = d.groupby('day').survival_percent.agg(n='count', mean='mean', sd='std').reset_index()
    d.to_csv(OUT / 'fig1f_survival_records.csv', index=False)
    summary.to_csv(OUT / 'fig1f_survival_summary.csv', index=False)
    plt.rcParams.update({'font.family': 'Times New Roman', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(5.06, 3.91))
    ax = fig.add_axes([.17, .145, .8, .78])
    for _, group in d.groupby('culture'):
        ax.plot(group.day, group.survival_percent, color='grey', alpha=.5, linewidth=1.6)
    ax.plot(summary.day, summary['mean'], color='forestgreen', linewidth=3)
    ax.set(xlim=(1, 20), ylim=(.0005, 3000), yscale='log')
    ax.set_xticks(range(2, 21, 2))
    ax.set_yticks([.001, .01, .1, 1, 10, 100, 1000], ['0.001', '0.01', '0.1', '1', '10', '100', '1000'])
    for name, spine in ax.spines.items():
        spine.set_visible(name in ('left', 'bottom'))
        spine.set_linewidth(1.1)
    ax.tick_params(which='major', labelsize=14, length=4, width=1.1, pad=2)
    ax.tick_params(which='minor', length=2, width=.6)
    ax.set_title('Evolution to Cefepime', fontsize=16, pad=5)
    ax.set_xlabel('Day', fontsize=16, labelpad=2)
    ax.set_ylabel('Survival (%)', fontsize=16, labelpad=2)
    target = OUT / 'trajectories'
    target.mkdir(exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(target / f'Fig1f_survival.{ext}', dpi=600)
    plt.close(fig)
    print(summary.tail(1).to_string(index=False))


if __name__ == '__main__':
    main()
