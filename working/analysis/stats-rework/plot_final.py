"""Dose-response, growth and survival figure panels.

Every tested IC50 panel: each culture as a point, the geometric mean as a line,
its 95% confidence interval as an error bar, and a bracket carrying the
Holm-adjusted p-value for each contrast that clears 0.05 (all adjusted p-values
are in out/final_statistics.csv and in the legend). Figure 5b-d and
Supplementary Figure 10a-c: the three technical series as points and the
geometric mean, no interval and no bracket (one culture per strain). Figure 5a:
six technical wells, mean and SD. Supplementary Figure 11: the recomputed
double-mutant MDK, all replicates shown, median and scaled MAD.

Styling follows the published panels: cohort colours and order, log axes on
whole decades, spine width 4, 15 x 10 inches, labelled and unlabelled
variants. Times New Roman stands in for Nimbus Roman. Output: PNG at 600 dpi
plus PDF and SVG with editable text, in out/figures-final/.
"""

import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
os.environ.setdefault("MPLCONFIGDIR", str(Path.home() / ".cache/matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy import stats as sps

from stats import split_strain
from fig5a_doubling import load_doubling_times
from manifest import panel_manifest, comparison_cohort
from report_results import p_text

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "figures-final"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Times New Roman"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

BRACKET_FIRST = 0.035  # axes fraction above the highest datum/interval
BRACKET_STEP = 0.075   # axes fraction; independent of the log data range
RENDERED = []

LAB = {"P": "MG", "PA": r"MG$^{\mathrm{AMI}}$", "PC": r"MG$^{\mathrm{CEF}}$",
       "PL": r"MG$^{\mathrm{LEV}}$", "PLA": r"MG$^{\mathrm{LEV,AMI}}$",
       "PLAC": r"MG$^{\mathrm{LEV,AMI,CEF}}$", "MG": "MG",
       "PCr": r"MG$^{\mathrm{CEF-R}}$",
       "ATEC": "Pb", "ATEC-C": r"Pb$^{\mathrm{CEF}}$", "ATEC-C-r": r"Pb$^{\mathrm{CEF-R}}$",
       "P-unt": r"MG$^{\mathrm{UNT}}$",
       "gata": r"MG$^{\mathrm{gatA}}$", "glvc": r"MG$^{\mathrm{glvC}}$",
       "hipa": r"MG$^{\mathrm{hipA}}$", "selb": r"MG$^{\mathrm{selB}}$",
       "rpoz": r"MG$^{\mathrm{rpoZ}}$", "ftsh": r"MG$^{\mathrm{ftsH}}$",
       "fime": r"MG$^{\mathrm{fimE}}$"}
COL = {"P": "grey", "PA": "blue", "PC": "forestgreen", "PL": "red",
       "PLA": "blue", "PLAC": "forestgreen", "MG": "grey",
       "PCr": "orchid", "ATEC": "grey", "ATEC-C": "forestgreen", "ATEC-C-r": "orchid",
       "P-unt": "firebrick",
       "gata": "red", "glvc": "orange", "hipa": "limegreen", "selb": "green",
       "rpoz": "blue", "ftsh": "purple", "fime": "brown"}
MUT = ["MG", "gata", "glvc", "hipa", "selb", "rpoz", "ftsh", "fime"]
DRUGS = ["Amikacin", "Levofloxacin", "Cefepime"]

# name, dataset, order, file stem, value, interval, brackets
PANELS = [
    ("Figure 2b-d", "fig2_paplpc", ["P", "PA", "PL", "PC"], "Fig2", "IC50", True, True),
    ("Figure 3b-d", "fig3_plac", ["P", "PL", "PLA", "PLAC"], "Fig3", "IC50", True, True),
    ("Supplementary Figure 3b-d", "supp3_pcr", ["P", "PCr"], "SuppFig3", "IC50", True, True),
    ("Supplementary Figure 4b-d", "supp4_atec", ["ATEC", "ATEC-C", "ATEC-C-r"], "SuppFig4",
     "IC50", True, True),
    ("Figure 5b-d", "fig5_mutants", MUT, "Fig5", "IC50", False, False),
    ("Supplementary Figure 10a-c", "fig5_mutants", MUT, "SuppFig10", "MIC", False, False),
    ("Figure 2b-d", "fig2_paplpc", ["P", "PA", "PL", "PC"], "SuppFig2", "MIC", True, False),
    ("Figure 3b-d", "fig3_plac", ["P", "PL", "PLA", "PLAC"], "SuppFig7", "MIC", True, False),
    ("Supplementary Figure 6a-c", "supp6_unt", ["P", "P-unt"], "SuppFig6", "IC50", True, True),
]


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=600 if ext == "png" else None,
                    bbox_inches="tight")
    RENDERED.append(stem)


def spread(n, x, s=0.36, gap=0.08):
    if n == 1:
        return np.array([float(x)])
    h = n // 2
    off = (np.concatenate([np.linspace(-s, -gap, h), np.linspace(gap, s, h)])
           if n % 2 == 0 else
           np.concatenate([np.linspace(-s, -gap, h), [0.0], np.linspace(gap, s, h)]))
    return x + off


def geo_ci(v):
    lg = np.log10(np.asarray(v, float))
    lg = lg[np.isfinite(lg)]
    n = len(lg)
    if n < 2:
        return (10 ** lg[0] if n else np.nan), np.nan, np.nan
    h = sps.t.ppf(0.975, n - 1) * lg.std(ddof=1) / np.sqrt(n)
    return 10 ** lg.mean(), 10 ** (lg.mean() - h), 10 ** (lg.mean() + h)


def nice_limits(values, n_brackets=0, ci_hi=None, ci_lo=None):
    """Whole-decade limits with headroom for the brackets."""
    v = np.asarray([x for x in values if np.isfinite(x) and x > 0], float)
    if not len(v):
        return 0.01, 1.0
    dmin, dmax = min(v.min(), ci_lo or v.min()), max(v.max(), ci_hi or 0)
    lo_e = int(np.floor(np.log10(dmin * 0.96)))
    reserved = BRACKET_FIRST + n_brackets * BRACKET_STEP + .025 if n_brackets else .06
    hi_e = int(np.ceil(lo_e + (np.log10(dmax) - lo_e + .04) / (1 - reserved)))
    hi_e = max(hi_e, lo_e + 1)
    return 10.0 ** lo_e, 10.0 ** hi_e


def style_log_axis(ax, ylim, present, labelled):
    ax.set_yscale("log")
    ax.set_ylim(ylim)
    ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10.0))
    ax.yaxis.set_minor_locator(mticker.LogLocator(
        base=10.0, subs=tuple(np.arange(2, 10) * 0.1) + tuple(range(2, 10))))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xlim(-0.7, len(present) - 0.3)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([LAB.get(g, g) for g in present] if labelled else [""] * len(present),
                       rotation=45, ha="center", fontsize=37)
    for l in ax.get_yticklabels():
        l.set_fontsize(37)
    for s in ax.spines.values():
        s.set_linewidth(4)
    ax.tick_params(axis="both", which="major", length=12, width=3)
    ax.tick_params(axis="both", which="minor", length=8, width=2)


def brackets_for(final, panel, drug, present):
    sel = final[(final.Panel == panel) & (final.value == "IC50") & (final.Antibiotic == drug)
                & (final.p_holm < 0.05)]
    rows = [r for r in sel.to_dict("records") if r["Group1"] in present and r["Group2"] in present]
    rows.sort(key=lambda r: abs(present.index(r["Group2"]) - present.index(r["Group1"])))
    return rows


def draw_panel(ax, sub, order, value, ylim, labelled, brackets, interval):
    present = [g for g in order if g in set(sub["group"])]
    tops = []
    for i, g in enumerate(present):
        v = sub[sub["group"] == g][value].dropna().values
        if not len(v):
            continue
        if interval and len(v) >= 10:
            ax.boxplot([v], positions=[i], widths=.48, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=COL.get(g, 'grey'), alpha=.12, linewidth=2),
                       medianprops=dict(color='grey', linewidth=2),
                       whiskerprops=dict(color='grey', linewidth=2),
                       capprops=dict(color='grey', linewidth=2), zorder=3, manage_ticks=False)
        ax.scatter(spread(len(v), i), v, s=170, facecolor=COL.get(g, "grey"),
                   edgecolor="black", linewidth=1.8, alpha=1, zorder=10)
        gm, lo, hi = geo_ci(v)
        ax.hlines(gm, i - 0.42, i + 0.42, color="black", linewidth=4, zorder=7)
        if interval and np.isfinite(lo) and lo > 0:
            ax.errorbar(i, gm, yerr=[[gm - lo], [hi - gm]], fmt="none", ecolor="black",
                        elinewidth=3, capsize=10, capthick=3, zorder=7)
            tops.append(hi)
        tops.append(v.max())
    style_log_axis(ax, ylim, present, labelled)
    loglo, loghi = np.log10(ylim)
    lvl = (np.log10(max(tops) if tops else ylim[0]) - loglo) / (loghi - loglo) + BRACKET_FIRST
    riser = 0.012
    for r in brackets:
        x1, x2 = sorted((present.index(r["Group1"]), present.index(r["Group2"])))
        if lvl + riser + .055 > 1:
            raise ValueError('Insufficient bracket headroom; no bracket may be omitted')
        ax.plot([x1, x1, x2, x2], [lvl, lvl+riser, lvl+riser, lvl],
                color="black", linewidth=2.5, zorder=12, transform=ax.get_xaxis_transform())
        p = r["p_holm"]
        txt = p_text(p)
        ax.text((x1 + x2) / 2, lvl + riser + .004, txt, ha="center", va="bottom",
                fontsize=36, transform=ax.get_xaxis_transform())
        lvl += BRACKET_STEP


def tested_panels(final):
    made = 0
    for panel, ds, order, stem, value, interval, use_brackets in PANELS:
        gf = comparison_cohort(ds, split_strain(pd.read_pickle(HERE / "cache" / f"{ds}.pkl")))
        gf["culturenumber"] = gf["culturenumber"].astype(str)
        for drug in DRUGS:
            sub = gf[gf["Antibiotic"] == drug]
            if sub.empty:
                continue
            present = [g for g in order if g in set(sub["group"])]
            br = brackets_for(final, panel, drug, present) if use_brackets else []
            vals = sub[sub["group"].isin(present)][value].dropna().values
            ci_hi = None
            ci_lo = None
            if interval:
                tops = [geo_ci(sub[sub["group"] == g][value].dropna().values)[2]
                        for g in present if (sub["group"] == g).sum() >= 2]
                ci_hi = max(t for t in tops if np.isfinite(t)) if tops else None
                lows = [geo_ci(sub[sub['group'] == g][value].dropna().values)[1]
                        for g in present if (sub['group'] == g).sum() >= 2]
                ci_lo = min(t for t in lows if np.isfinite(t)) if lows else None
            ylim = nice_limits(vals, len(br), ci_hi, ci_lo)
            for labelled in (True, False):
                fig, ax = plt.subplots(figsize=(15, 10))
                draw_panel(ax, sub, order, value, ylim, labelled, br, interval)
                sym = "IC$_{50}$" if value == "IC50" else "MIC"
                ax.set_ylabel(f"{drug}\n{sym} (μg/mL)", fontsize=40)
                fig.subplots_adjust(left=0.2, right=0.9, top=0.9)
                save(fig, f"{stem}_{drug}_{value}{'_labels' if labelled else ''}")
                plt.close(fig)
                made += 1
    return made


def figure5a():
    g = load_doubling_times()
    order = ["MG", "gatA", "glvC", "hipA", "selB", "rpoZ", "ftsH", "fimE"]
    key = {"MG": "MG", "gatA": "gata", "glvC": "glvc", "hipA": "hipa", "selB": "selb",
           "rpoZ": "rpoz", "ftsH": "ftsh", "fimE": "fime"}
    for labelled in (True, False):
        fig, ax = plt.subplots(figsize=(15, 10))
        for i, s in enumerate(order):
            v = g[g["Culture"] == s]["doubling_time_min"].dropna().values
            ax.scatter(spread(len(v), i), v, s=170, facecolor="grey", edgecolor="black",
                       linewidth=1.8, alpha=1, zorder=10)
            ax.hlines(v.mean(), i - 0.42, i + 0.42, color="black", linewidth=4, zorder=7)
            ax.errorbar(i, v.mean(), yerr=v.std(ddof=1), fmt="none", ecolor="black",
                        elinewidth=3, capsize=10, capthick=3, zorder=7)
        ax.set_ylim(0, 25)
        ax.set_xlim(-0.7, len(order) - 0.3)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([LAB[key[s]] for s in order] if labelled else [""] * len(order),
                           rotation=45, ha="center", fontsize=37)
        for l in ax.get_yticklabels():
            l.set_fontsize(37)
        for sp in ax.spines.values():
            sp.set_linewidth(4)
        ax.tick_params(axis="both", which="major", length=12, width=3)
        ax.set_ylabel("Doubling time (min)", fontsize=40)
        fig.subplots_adjust(left=0.2, right=0.9, top=0.9)
        save(fig, f"Fig5a_doubling_time{'_labels' if labelled else ''}")
        plt.close(fig)
    return 2


def draw_mdk(rep, summ, order, col, lab, stem, conditional=False):
    fig, ax = plt.subplots(figsize=(10, 12))
    offsets = dict(zip(order, np.linspace(-.25, .25, len(order))))
    omitted = 0
    for name in order:
        s = summ[summ.Name == name].sort_values('Time')
        med, mad = s['median'].values, s.mad_scaled.values
        t = s.Time.values + offsets[name]
        lower_ok = med - mad > 0
        omitted += int((~lower_ok).sum())
        ax.errorbar(t, med, yerr=[np.where(lower_ok, mad, 0), mad], color=col[name],
                    linewidth=3, marker='o', markersize=7, capsize=5, capthick=2,
                    elinewidth=2, alpha=.8, label=lab[name], zorder=4)
        r = rep[rep.Name == name].copy()
        x = np.empty(len(r))
        for time in r.Time.unique():
            mask = (r.Time == time).values
            x[mask] = spread(mask.sum(), time + offsets[name], s=.035, gap=.012)
        cen = r.below_detection.to_numpy() if 'below_detection' in r else np.zeros(len(r), bool)
        values = r.fraction_plot.to_numpy() if 'fraction_plot' in r else r.fraction_observed.to_numpy()
        ax.scatter(x[~cen], values[~cen], s=55, facecolor=col[name], edgecolor='black',
                   linewidth=.8, alpha=1, zorder=8)
        if cen.any():
            ax.scatter(x[cen], values[cen], s=140, facecolor='white', edgecolor=col[name],
                       linewidth=2, marker='v', zorder=10)
        if 'mad_scaled_max' in s:
            # Dashed extension shows the additional spread allowed by the censored observation.
            for x0, mid, low_mad, high_mad in zip(t, med, s.mad_scaled_min, s.mad_scaled_max):
                if high_mad > low_mad:
                    ax.plot([x0, x0], [mid + low_mad, mid + high_mad], ls='--', color=col[name], lw=2)
                    if mid - high_mad > 0:
                        ax.plot([x0, x0], [mid - high_mad, mid - low_mad], ls='--', color=col[name], lw=2)
    ax.set_yscale('log')
    ax.set_ylim(1e-9, 1.5)
    ax.set_xlim(-.5, 7.5)
    ax.set_xticks([0, 1, 2, 3, 5, 7])
    ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10.0))
    ax.set_xlabel('Time (h)', fontsize=36)
    ax.set_ylabel('Survivor fraction', fontsize=36)
    for label in ax.get_yticklabels() + ax.get_xticklabels():
        label.set_fontsize(30)
    for sp in ax.spines.values():
        sp.set_linewidth(4)
    ax.tick_params(axis='both', which='major', length=12, width=3)
    ax.legend(fontsize=24, frameon=False, loc='lower center', bbox_to_anchor=(.5, 1.02), ncol=2)
    notes = ['Horizontal offsets separate strains and replicates.']
    if conditional:
        notes += ['10 µL plated at every measured time.',
                  'Open triangle: below detection. Dashed bar: allowed MAD range.']
    if omitted:
        notes += ['Lower MAD arms reaching zero are omitted on the log axis.']
    fig.text(.18, .015, '\n'.join(notes), fontsize=16, ha='left', va='bottom')
    fig.subplots_adjust(left=.2, right=.95, top=.82, bottom=.20 if conditional else .14)
    save(fig, stem)
    plt.close(fig)


def mdk_panels():
    rep = pd.read_csv(HERE / 'out/supp11_recomputed.csv')
    summ = pd.read_csv(HERE / 'out/supp11_summary.csv')
    order = ['hipA', 'hipA_selB', 'hipA_rpoZ', 'hipA_ftsH', 'wt']
    col = {'wt': 'black', 'hipA_rpoZ': 'blueviolet', 'hipA_ftsH': 'steelblue',
           'hipA': 'violet', 'hipA_selB': 'lightsalmon'}
    lab = {'wt': 'MG', 'hipA_rpoZ': r'MG$^{\mathrm{hipA+rpoZ}}$',
           'hipA_ftsH': r'MG$^{\mathrm{hipA+ftsH}}$', 'hipA': r'MG$^{\mathrm{hipA}}$',
           'hipA_selB': r'MG$^{\mathrm{hipA+selB}}$'}
    draw_mdk(rep, summ, order, col, lab, 'SuppFig11_MDK_double_mutants', conditional=True)
    for stem, order, col, lab in [
        ('Fig5e_MDK', ['gata', 'glvc', 'hipa', 'selb', 'rpoz', 'ftsh', 'fimE', 'wt'],
         dict(zip(['gata','glvc','hipa','selb','rpoz','ftsh','fimE','wt'],
                  ['teal','blue','violet','maroon','orange','purple','red','black'])),
         {**{g: LAB[g] for g in ['gata','glvc','hipa','selb','rpoz','ftsh']}, 'fimE': LAB['fime'], 'wt': 'MG'}),
        ('SuppFig10d_MDK', ['hipa', 'plac', 'prs', 'pc', 'wt'],
         {'hipa':'violet','plac':'purple','prs':'teal','pc':'blue','wt':'black'},
         {'hipa': LAB['hipa'], 'plac': LAB['PLAC']+'-7', 'prs':r'MG$^{\mathrm{prs}}$',
          'pc':LAB['PC']+'-5', 'wt':'MG'}),
    ]:
        rep = pd.read_csv(HERE / 'out' / f'{stem}_records.csv')
        summ = pd.read_csv(HERE / 'out' / f'{stem}_summary.csv')
        draw_mdk(rep, summ, order, col, lab, stem)
    return 3


def main():
    final = pd.read_csv(HERE / "out" / "final_statistics.csv")
    n = tested_panels(final) + figure5a() + mdk_panels()
    manifest = panel_manifest()
    expected = {Path(x).name for x in manifest[manifest.status == 'standalone plot'].output}
    actual = {x.removesuffix('_labels') for x in RENDERED}
    if expected != actual:
        raise ValueError(f'Render coverage mismatch: missing {expected-actual}, extra {actual-expected}')
    pd.DataFrame({'stem': RENDERED}).to_csv(HERE / 'out/render_manifest.csv', index=False)
    print(f"wrote {n} panels (x3 formats) to {OUT}")


if __name__ == "__main__":
    main()
