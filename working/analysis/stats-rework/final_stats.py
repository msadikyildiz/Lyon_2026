"""Statistical analysis of dose-response, growth and survival measurements.

The 17 August comparison table is retained under superseded/out-17aug.
Current conventions and figure coverage are described in the repository README.

Experimental design:
  * P1-P10 are the day-0 ancestors of every same-numbered culture in every arm,
    so contrasts among P, PA, PC, PL, PLA, PLAC are within-lineage and PAIRED,
    aligned on culture number. (The published code paired by row position.)
  * Figure 5b-d and Supplementary Figure 10 triplicates are three technical
    dose-response series from one overnight culture per strain: biological
    n = 1, so no inferential test. Descriptive tables only.
  * The other panels have no established lineage link between groups and use
    Welch's unpaired t-test.

Estimand: log10 IC50 (primary). Effects are ratios of geometric means with
pointwise 95% CIs. MIC is descriptive (geometric mean and CI, no test).
Correction: Holm within panel and drug. Sensitivity:
BH and Bonferroni per drug, Holm per figure, the other test (Welch or paired),
raw scale, exact sign-flip or permutation, nested bootstrap over cultures and
curve fits, and Figure 3 with lineages 2, 3, 8 excluded.

Run:  python final_stats.py
Reads cache/*.pkl (built by build_all.py). Writes out/*.csv.
"""

import itertools
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from stats import split_strain
from fig5a_doubling import load_doubling_times
from manifest import EXPECTED_IDS

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
SEED = 20260915
RNG = np.random.default_rng(SEED)
NBOOT = 4000
MC_PERM = 100_000
SUSPECT_LINEAGES = {"2", "3", "8"}   # mutL hypermutators, possible contamination

# panel, dataset, contrasts, design, also shown in
PANELS = [
    ("Figure 2b-d", "fig2_paplpc", [("P", "PA"), ("P", "PC"), ("P", "PL")],
     "paired", "Supplementary Figure 2"),
    ("Figure 3b-d", "fig3_plac",
     list(itertools.combinations(["P", "PL", "PLA", "PLAC"], 2)),
     "paired", "Supplementary Figure 7"),
    ("Supplementary Figure 3b-d", "supp3_pcr", [("P", "PCr")], "welch", ""),
    ("Supplementary Figure 4b-d", "supp4_atec", [("ATEC", "ATEC-C")], "welch", ""),
    ("Supplementary Figure 6a-c", "supp6_unt", [("P", "P-unt")], "welch", ""),
]
DRUG_ABBR = {"Levofloxacin": "LEV", "Amikacin": "AMI", "Cefepime": "CEF"}
MUTANTS = ["gata", "glvc", "hipa", "selb", "rpoz", "ftsh", "fime"]
MUT_LABEL = {"gata": "gatA", "glvc": "glvC", "hipa": "hipA", "selb": "selB",
             "rpoz": "rpoZ", "ftsh": "ftsH", "fime": "fimE", "MG": "MG"}


# --------------------------------------------------------------------------- io
def load(ds):
    gf = split_strain(pd.read_pickle(HERE / "cache" / f"{ds}.pkl"))
    gf["culturenumber"] = gf["culturenumber"].astype(str)
    return gf


def draws(cell):
    """The 1000 bootstrap IC50 fits stored per curve, as log10."""
    if isinstance(cell, str):
        cell = json.loads(cell.replace("nan", "NaN"))
    v = np.asarray(cell, float)
    if v.ndim != 1 or not len(v) or not np.all(np.isfinite(v) & (v > 0)):
        raise ValueError('Fit draws must be a nonempty finite positive vector')
    return np.log10(v)


def require_positive(values, label='measurements'):
    v = np.asarray(values, float)
    if not len(v) or not np.all(np.isfinite(v) & (v > 0)):
        raise ValueError(f'{label}: expected finite positive measurements')


def aligned(a, b):
    """Reject lost or duplicated pairs; never take an implicit intersection."""
    if not a.index.is_unique or not b.index.is_unique:
        raise ValueError('Duplicate culture IDs')
    if set(a.index) != set(b.index):
        raise ValueError(f'Incomplete pairs: {sorted(set(a.index) ^ set(b.index))}')
    require_positive(a); require_positive(b)
    keys = sorted(a.index, key=int)
    return keys, a.reindex(keys).astype(float).values, b.reindex(keys).astype(float).values


def validate_inputs(datasets):
    checks = []
    for ds, gf in datasets.items():
        expected = EXPECTED_IDS[ds]
        for drug in DRUG_ABBR:
            d = gf[gf.Antibiotic == drug]
            if set(d['group']) != set(expected):
                raise ValueError(f'{ds}/{drug}: unexpected or missing groups')
            for group, ids in expected.items():
                g = d[d['group'] == group]
                if not g.culturenumber.is_unique or set(g.culturenumber) != set(ids):
                    raise ValueError(f'{ds}/{drug}/{group}: unexpected, missing, or duplicate culture IDs')
                checks.append(('expected unique culture IDs', ds, drug, group, True))
                for value in ['IC50', 'MIC']:
                    require_positive(g[value], f'{ds}/{drug}/{group}/{value}')
                    checks.append((f'finite positive {value}', ds, drug, group, True))
                if ds != 'fig5_mutants':
                    for cell in g.ic50_bootstrap:
                        if len(draws(cell)) != 1000:
                            raise ValueError(f'{ds}/{drug}/{group}: expected 1000 stored fit draws')
                    checks.append(('1000 valid fit draws per culture', ds, drug, group, True))
    return checks


# ------------------------------------------------------------------ contrasts
def geo_summary(values):
    lg = np.log10(np.asarray(values, float))
    lg = lg[np.isfinite(lg)]
    n = len(lg)
    if n == 0:
        return dict(n=0, geomean=np.nan, lo=np.nan, hi=np.nan, min=np.nan, max=np.nan)
    if n == 1:
        return dict(n=1, geomean=10 ** lg[0], lo=np.nan, hi=np.nan,
                    min=10 ** lg[0], max=10 ** lg[0])
    h = stats.t.ppf(0.975, n - 1) * lg.std(ddof=1) / np.sqrt(n)
    return dict(n=n, geomean=10 ** lg.mean(), lo=10 ** (lg.mean() - h),
                hi=10 ** (lg.mean() + h), min=10 ** lg.min(), max=10 ** lg.max())


def paired(a, b):
    """a, b: Series indexed by culture number. Returns log10 paired contrast b/a."""
    common, av, bv = aligned(a, b)
    av, bv = np.log10(av), np.log10(bv)
    d = bv - av
    n = len(d)
    if n < 2:
        return dict(n_pairs=n, ratio=np.nan, lo=np.nan, hi=np.nan, p=np.nan, t=np.nan, df=np.nan)
    t, p = stats.ttest_rel(bv, av)
    half = stats.t.ppf(0.975, n - 1) * d.std(ddof=1) / np.sqrt(n)
    return dict(n_pairs=n, ratio=10 ** d.mean(), lo=10 ** (d.mean() - half),
                hi=10 ** (d.mean() + half), p=p, t=t, df=n - 1,
                r=(np.corrcoef(av, bv)[0, 1] if n > 2 else np.nan), _d=d,
                _av=av, _bv=bv, _keys=common)


def welch(a, b):
    av = np.log10(a.dropna().astype(float).values)
    bv = np.log10(b.dropna().astype(float).values)
    av, bv = av[np.isfinite(av)], bv[np.isfinite(bv)]
    if len(av) < 2 or len(bv) < 2:
        return dict(n_pairs=np.nan, ratio=np.nan, lo=np.nan, hi=np.nan, p=np.nan, t=np.nan, df=np.nan)
    t, p = stats.ttest_ind(bv, av, equal_var=False)
    m = bv.mean() - av.mean()
    va, vb = av.var(ddof=1) / len(av), bv.var(ddof=1) / len(bv)
    dfw = (va + vb) ** 2 / (va ** 2 / (len(av) - 1) + vb ** 2 / (len(bv) - 1))
    half = stats.t.ppf(0.975, dfw) * np.sqrt(va + vb)
    return dict(n_pairs=np.nan, ratio=10 ** m, lo=10 ** (m - half), hi=10 ** (m + half),
                p=p, t=t, df=dfw, _av=av, _bv=bv)


def raw_p(a, b, design):
    """The same test on the raw (unlogged) scale, for the sensitivity table."""
    if design == "paired":
        _, av, bv = aligned(a, b)
        return stats.ttest_rel(bv, av).pvalue if len(av) > 1 else np.nan
    av, bv = a.dropna().astype(float).values, b.dropna().astype(float).values
    return stats.ttest_ind(bv, av, equal_var=False).pvalue if min(len(av), len(bv)) > 1 else np.nan


def signflip_p(d):
    """Exact two-sided sign-flip test on paired log differences."""
    d = np.asarray(d, float)
    n = len(d)
    obs = abs(d.mean())
    signs = np.array(list(itertools.product([-1, 1], repeat=n)))
    means = np.abs((signs * d).mean(axis=1))
    return float((means >= obs - 1e-12).mean())


def perm_p(av, bv):
    """Exact (or Monte Carlo) two-sided permutation test on the difference of means."""
    pooled = np.concatenate([av, bv])
    na, n = len(av), len(pooled)
    obs = abs(bv.mean() - av.mean())
    from math import comb
    if comb(n, na) <= 200_000:
        cnt = tot = 0
        for idx in itertools.combinations(range(n), na):
            mask = np.zeros(n, bool); mask[list(idx)] = True
            diff = abs(pooled[~mask].mean() - pooled[mask].mean())
            cnt += diff >= obs - 1e-12; tot += 1
        return cnt / tot
    cnt = 0
    for _ in range(MC_PERM):
        perm = RNG.permutation(pooled)
        cnt += abs(perm[na:].mean() - perm[:na].mean()) >= obs - 1e-12
    return (cnt + 1) / (MC_PERM + 1)


def nested_boot(da, db, design):
    """Return a fit-perturbation percentile range, without a bootstrap test.

    da, db: lists of log10 draw arrays, one per culture; for 'paired' they are aligned.
    """
    na, nb = len(da), len(db)
    if not na or not nb or (design == 'paired' and na != nb):
        raise ValueError('Missing culture arrays or unequal paired array counts')
    if any(np.asarray(x).ndim != 1 or not len(x) or not np.isfinite(x).all() for x in [*da, *db]):
        raise ValueError('Empty or invalid fit-draw array; pairs must remain intact')
    diffs = np.empty(NBOOT)
    for i in range(NBOOT):
        if design == "paired":
            idx = RNG.integers(0, na, na)
            xa = np.array([RNG.choice(da[j]) for j in idx])
            xb = np.array([RNG.choice(db[j]) for j in idx])
        else:
            ia = RNG.integers(0, na, na); ib = RNG.integers(0, nb, nb)
            xa = np.array([RNG.choice(da[j]) for j in ia])
            xb = np.array([RNG.choice(db[j]) for j in ib])
        diffs[i] = xb.mean() - xa.mean()
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return 10 ** lo, 10 ** hi


# ------------------------------------------------------------------ pipeline
def contrast_rows(datasets=None):
    rows, groups = [], []
    checks = []
    for panel, ds, comps, design, also in PANELS:
        gf = datasets[ds] if datasets is not None else load(ds)
        present = set(gf["group"])
        for drug in ["Levofloxacin", "Amikacin", "Cefepime"]:
            d = gf[gf["Antibiotic"] == drug]
            # group summaries, both measures
            for g in sorted(present):
                dg = d[d["group"] == g]
                if not len(dg):
                    continue
                # ID uniqueness
                dup = dg["culturenumber"].duplicated().any()
                checks.append(("unique culture IDs", panel, drug, g, not dup))
                for value in ("IC50", "MIC"):
                    s = geo_summary(dg[value])
                    groups.append(dict(Panel=panel, Antibiotic=drug, Group=g, value=value,
                                       cultures=",".join(dg["culturenumber"]),
                                       values=";".join(repr(float(x)) for x in dg[value]),
                                       **s))
            for g1, g2 in comps:
                if g1 not in present or g2 not in present:
                    continue
                a = d[d["group"] == g1].set_index("culturenumber")
                b = d[d["group"] == g2].set_index("culturenumber")
                for value in ("IC50", "MIC"):
                    av, bv = a[value], b[value]
                    res = paired(av, bv) if design == "paired" else welch(av, bv)
                    row = dict(Panel=panel, also_shown_in=also, value=value, Antibiotic=drug,
                               contrast_id=f"{DRUG_ABBR[drug]}:{g1}-vs-{g2}",
                               Group1=g1, Group2=g2,
                               n1=int(av.notna().sum()), n2=int(bv.notna().sum()),
                               n_pairs=res.get("n_pairs"),
                               geomean1=geo_summary(av)["geomean"], geomean2=geo_summary(bv)["geomean"],
                               ratio=res["ratio"], ratio_lo=res["lo"], ratio_hi=res["hi"],
                               test=(("paired t on log10, aligned on culture"
                                      if design == "paired" else "Welch t on log10")
                                     if value == 'IC50' else 'descriptive'),
                               interval_method='paired log-ratio t interval' if design == 'paired'
                                               else 'Welch log-ratio t interval',
                               t=res.get("t") if value == 'IC50' else np.nan,
                               df=res.get("df") if value == 'IC50' else np.nan,
                               r_pairs=res.get("r", np.nan) if value == 'IC50' else np.nan)
                    if design == 'paired':
                        checks.append((f'complete {value} pairs', panel, drug, f'{g1}-{g2}',
                                       res['n_pairs'] == len(av) == len(bv)))
                    if value == "IC50":
                        row["p"] = res["p"]
                        # sensitivities
                        common = sorted(set(av.index) & set(bv.index), key=int)
                        other = (welch(av, bv) if design == "paired"
                                 else paired(av.loc[common], bv.loc[common]))
                        row["p_other_test"] = other["p"]
                        row["other_test"] = ('Welch' if design == 'paired' else
                                             'paired on shared numeric IDs; matching unconfirmed')
                        row['other_test_ids'] = ','.join(common) if design != 'paired' else ''
                        row["p_raw_scale"] = raw_p(av, bv, design)
                        if design == "paired" and res["n_pairs"] >= 2:
                            row["p_exact"] = signflip_p(res["_d"])
                            row["exact_test"] = "sign-flip"
                            da = [draws(a.loc[k, "ic50_bootstrap"]) for k in res["_keys"]]
                            db = [draws(b.loc[k, "ic50_bootstrap"]) for k in res["_keys"]]
                        elif design != "paired" and not np.isnan(res["p"]):
                            row["p_exact"] = perm_p(res["_av"], res["_bv"])
                            row["exact_test"] = "permutation"
                            da = [draws(x) for x in a["ic50_bootstrap"]]
                            db = [draws(x) for x in b["ic50_bootstrap"]]
                        else:
                            da = db = None
                        if da:
                            lo, hi = nested_boot(da, db, design)
                            row.update(fit_range_lo=lo, fit_range_hi=hi,
                                       fit_range_method='culture resampling and fit perturbation; empirical 2.5-97.5 percentiles',
                                       fit_resamples=NBOOT, fit_seed=SEED)
                        if panel == "Figure 3b-d":
                            keep = [k for k in a.index if k not in SUSPECT_LINEAGES]
                            r2 = paired(av.loc[[k for k in keep if k in av.index]],
                                        bv.loc[[k for k in keep if k in bv.index]])
                            row["p_excl_lineages_238"] = r2["p"]
                            row["ratio_excl_lineages_238"] = r2["ratio"]
                    rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(groups), checks


def correct(t):
    """Holm per panel x drug (primary), BH and Bonferroni per panel x drug, Holm per panel."""
    ic = t["value"] == "IC50"
    for col in ("p_holm", "p_bh", "p_bonferroni", "p_holm_figure", "p_holm_excl_238"):
        t[col] = np.nan
    t["k"] = np.nan
    for (_, _), idx in t[ic].groupby(["Panel", "Antibiotic"]).groups.items():
        sub = t.loc[idx, "p"]; ok = sub.notna()
        t.loc[idx, "k"] = int(ok.sum())
        if ok.sum():
            t.loc[sub[ok].index, "p_holm"] = multipletests(sub[ok], method="holm")[1]
            t.loc[sub[ok].index, "p_bh"] = multipletests(sub[ok], method="fdr_bh")[1]
            t.loc[sub[ok].index, "p_bonferroni"] = np.minimum(sub[ok] * ok.sum(), 1)
        if "p_excl_lineages_238" in t:
            s2 = t.loc[idx, "p_excl_lineages_238"]; ok2 = s2.notna()
            if ok2.sum():
                t.loc[s2[ok2].index, "p_holm_excl_238"] = multipletests(s2[ok2], method="holm")[1]
    for _, idx in t[ic].groupby("Panel").groups.items():
        sub = t.loc[idx, "p"]; ok = sub.notna()
        if ok.sum():
            t.loc[sub[ok].index, "p_holm_figure"] = multipletests(sub[ok], method="holm")[1]
    for raw_col in ['p_other_test', 'p_raw_scale', 'p_exact']:
        adjusted_col = raw_col + '_holm'
        t[adjusted_col] = np.nan
        for _, idx in t[ic].groupby(['Panel', 'Antibiotic']).groups.items():
            valid = t.loc[idx, raw_col].dropna()
            if len(valid):
                t.loc[valid.index, adjusted_col] = multipletests(valid, method='holm')[1]
    t['sensitivity_family'] = np.where(ic, 'same panel and antibiotic as primary', '')
    t['significant'] = pd.Series(pd.NA, index=t.index, dtype='boolean')
    t.loc[ic, 'significant'] = t.loc[ic, 'p_holm'] < 0.05
    return t


def figure5_descriptive():
    gf = load("fig5_mutants")
    rows = []
    for value in ("IC50", "MIC"):
        for drug in ["Levofloxacin", "Amikacin", "Cefepime"]:
            d = gf[gf["Antibiotic"] == drug]
            mg = geo_summary(d[d["group"] == "MG"][value])
            for g in ["MG"] + MUTANTS:
                dg = d[d["group"] == g][value].astype(float)
                s = geo_summary(dg)
                rows.append(dict(Panel="Figure 5b-d / Supplementary Figure 10a-c", value=value,
                                 Antibiotic=drug, Strain=MUT_LABEL[g],
                                 replicates="three technical dose-response series, one overnight culture",
                                 values=";".join(repr(float(x)) for x in dg),
                                 n_technical=s["n"], geomean=s["geomean"], min=s["min"], max=s["max"],
                                 ratio_vs_MG=s["geomean"] / mg["geomean"],
                                 ratio_min=s["min"] / mg["geomean"], ratio_max=s["max"] / mg["geomean"]))
    return pd.DataFrame(rows)


def figure5a_descriptive():
    g = load_doubling_times()
    rows = []
    mg = g[g["Culture"] == "MG"]["doubling_time_min"].dropna()
    for strain in ["MG", "gatA", "glvC", "selB", "hipA", "rpoZ", "ftsH", "fimE"]:
        v = g[g["Culture"] == strain]["doubling_time_min"].dropna()
        rows.append(dict(Panel="Figure 5a", Strain=strain,
                         replicates="six technical wells, one overnight culture",
                         values=";".join(repr(float(x)) for x in v), n_technical=len(v),
                         mean_min=v.mean(), sd_min=v.std(ddof=1), min=v.min(), max=v.max(),
                         ratio_vs_MG=v.mean() / mg.mean()))
    return pd.DataFrame(rows)


def what_changed(t):
    """Published call (uncorrected, raw scale, positional pairing) versus the primary call."""
    from validate import CHECKS
    pub = {}
    for ds, mode, value, label, cases in CHECKS:
        for drug, g1, g2, p in cases:
            pub[(ds, drug, g1, g2)] = p
    ds_of = {p[0]: p[1] for p in PANELS}
    rows = []
    for _, r in t[t.value == "IC50"].iterrows():
        key = (ds_of[r.Panel], r.Antibiotic, r.Group1, r.Group2)
        p_pub = pub.get(key, np.nan)
        rows.append(dict(Panel=r.Panel, contrast_id=r.contrast_id, ratio=r.ratio,
                         p_published=p_pub, published_call=(p_pub < 0.05) if pd.notna(p_pub) else None,
                         p=r.p, p_holm=r.p_holm, revised_call=bool(r.significant),
                         changed=(pd.notna(p_pub) and (p_pub < 0.05) != bool(r.significant))))
    return pd.DataFrame(rows)


def versus_17aug(t):
    old = pd.read_csv(HERE / "superseded" / "out-17aug" / "final_statistics_17aug.csv")
    old = old[old.value == "IC50"]
    old['Panel'] = old.Panel.replace({'Supplementary Figure 5a / 6a': 'Supplementary Figure 6a-c'})
    used = set()
    rows = []
    for _, r in t[t.value == "IC50"].iterrows():
        o = old[(old.Panel == r.Panel) & (old.Antibiotic == r.Antibiotic)
                & (old.Group1 == r.Group1) & (old.Group2 == r.Group2)]
        if len(o) != 1:
            raise ValueError(f'Historical comparison missing or duplicated: {r.Panel}/{r.contrast_id}')
        used.add(o.index[0])
        o = o.iloc[0]
        rows.append(dict(Panel=r.Panel, contrast_id=r.contrast_id,
                         disposition='matched',
                         test_17aug=o.test, q_bh_17aug=o.q_bh, call_17aug=bool(o.sig_bh),
                         test_now=r.test, p_holm_now=r.p_holm, call_now=bool(r.significant),
                         changed=bool(o.sig_bh) != bool(r.significant)))
    retired = old.loc[~old.index.isin(used)]
    if set(retired.Panel) - {'Figure 5b-d'}:
        raise ValueError('Unexpected unmatched historical comparisons')
    for _, o in retired.iterrows():
        rows.append(dict(Panel=o.Panel, contrast_id=f'{DRUG_ABBR[o.Antibiotic]}:{o.Group1}-vs-{o.Group2}',
                         disposition='retired: biological n=1, descriptive reporting',
                         test_17aug=o.test, q_bh_17aug=o.q_bh, call_17aug=bool(o.sig_bh),
                         test_now='descriptive', p_holm_now=np.nan, call_now=None, changed=False))
    return pd.DataFrame(rows)


def exclusion_checks():
    """Declared exclusions are present in the cached analysis tables."""
    out = []
    d3 = pd.read_pickle(HERE / "cache" / "fig3_plac__df_analysis.pkl")
    d3["Experiment"] = d3["Experiment"].astype(str)
    out.append(("PLAC7/PLAC10 absent from experiment 1", "Figure 3b-d", "", "",
                not ((d3.Experiment == "1") & d3.Strain.isin(["PLAC7", "PLAC10"])).any()))
    d4 = pd.read_pickle(HERE / "cache" / "supp4_atec__df_analysis.pkl")
    d4["Experiment"] = d4["Experiment"].astype(str)
    bad = ((d4.Experiment == "2") & (d4.Antibiotic == "Cefepime")
           & d4.Strain.isin(["ATEC1", "ATEC2", "ATEC4", "ATEC-C3", "ATEC-C-r1"])).any()
    out.append(("experiment-2 cefepime failures excluded (incl. ATEC-C-r1)",
                "Supplementary Figure 4b-d", "", "", not bad))
    keep = ((d4.Experiment == "2") & (d4.Antibiotic != "Cefepime")
            & d4.Strain.isin(["ATEC1", "ATEC2", "ATEC4", "ATEC-C3"])).sum() > 0
    out.append(("experiment-2 levofloxacin/amikacin data retained for those cultures",
                "Supplementary Figure 4b-d", "", "", keep))
    return out


def main():
    global RNG
    RNG = np.random.default_rng(SEED)
    datasets = {ds: load(ds) for ds in EXPECTED_IDS}
    checks = validate_inputs(datasets) + exclusion_checks()
    if any(not c[-1] for c in checks):
        raise ValueError(f'Input validation failed: {[c for c in checks if not c[-1]]}')
    t, groups, pair_checks = contrast_rows(datasets)
    t = correct(t)
    checks += pair_checks
    ic = t[t.value == 'IC50']
    checks += [('36 primary IC50 contrasts', '', '', '', len(ic) == 36),
               ('all primary results finite', '', '', '',
                np.isfinite(ic[['ratio', 'ratio_lo', 'ratio_hi', 'p', 'p_holm']]).all().all())]
    f5, f5a = figure5_descriptive(), figure5a_descriptive()
    wc, v17 = what_changed(t), versus_17aug(t)
    checks.append(('36 matched August comparisons', '', '', '', (v17.disposition == 'matched').sum() == 36))
    if any(not c[-1] for c in checks):
        raise ValueError(f'Output validation failed: {[c for c in checks if not c[-1]]}')
    primary_cols = ['Panel', 'also_shown_in', 'value', 'Antibiotic', 'contrast_id', 'Group1', 'Group2',
                    'n1', 'n2', 'n_pairs', 'geomean1', 'geomean2', 'ratio', 'ratio_lo', 'ratio_hi',
                    'test', 'interval_method', 'p', 'k', 'p_holm', 'significant', 'p_bh',
                    'p_bonferroni', 'p_holm_figure']
    outputs = {'final_statistics.csv': t[primary_cols], 'sensitivity_all.csv': t,
               'group_summaries.csv': groups, 'fig5_descriptive.csv': f5,
               'fig5a_descriptive.csv': f5a, 'fig5a_doubling_times_per_well.csv': load_doubling_times(),
               'what_changes_vs_published.csv': wc, 'what_changes_vs_17aug_report.csv': v17,
               'acceptance_checks.csv': pd.DataFrame(checks, columns=['check', 'panel', 'drug', 'item', 'passed'])}
    # Serialize all tables only after validation, then replace the existing files.
    from tempfile import TemporaryDirectory
    OUT.mkdir(exist_ok=True)
    with TemporaryDirectory(dir=OUT.parent) as temp:
        staging = Path(temp)
        for name, table in outputs.items():
            table.to_csv(staging / name, index=False)
        for name in outputs:
            (staging / name).replace(OUT / name)
    print(f'{len(checks)}/{len(checks)} acceptance checks passed; {len(ic)} primary IC50 contrasts')
    print(ic[['Panel', 'Antibiotic', 'Group1', 'Group2', 'ratio', 'p_holm']].to_string(index=False))
    return outputs


if __name__ == '__main__':
    main()
