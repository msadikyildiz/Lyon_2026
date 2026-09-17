"""Shared helpers and superseded August statistics for historical reproduction.

The current primary analysis is final_stats.py: Figures 2/3 use culture-ID-paired
log10 IC50 tests; Supplementary Figures 3/4/6 use Welch tests; Holm adjustment is
by panel and drug. See RESULTS_v2.md and the repository README.

compare() retains the old raw-scale modes for validate.py. panel_tests() and
its design_mode() reproduce the superseded August exploratory analysis only;
they do not select the current study design. The August all-Welch decision
predated the incorporation of Adam's culture-matching clarification and is
superseded. Its raw arithmetic means and exploratory corrections are historical.
"""

import itertools
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def split_strain(df):
    """Add the notebook's ``group`` and ``culturenumber`` columns."""
    out = df.copy()
    out["group"] = out["Strain"].apply(lambda x: re.split(r"[1-9]", str(x))[0])
    out["culturenumber"] = out["Strain"].apply(
        lambda x: (m.group() if (m := re.search(r"(\d+)", str(x))) else None))
    return out


def ci95(values):
    """Mean and 95% confidence interval half-width of the mean."""
    v = pd.Series(values).dropna().astype(float)
    n = len(v)
    if n < 2:
        return (v.mean() if n else np.nan), np.nan, n
    half = stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n)
    return v.mean(), half, n


def compare(a, b, mode, keys_a=None, keys_b=None):
    """Return (t, p, n_a, n_b, n_pairs) for one comparison.

    mode is 'welch', 'paired_fixed' (align on culture number) or
    'paired_published' (reproduce the positional-pairing bug).
    """
    a = pd.Series(list(a), index=None if keys_a is None else list(keys_a))
    b = pd.Series(list(b), index=None if keys_b is None else list(keys_b))
    if mode == "welch":
        av, bv = a.dropna(), b.dropna()
        if len(av) < 2 or len(bv) < 2:
            return np.nan, np.nan, len(av), len(bv), np.nan
        t, p = stats.ttest_ind(av, bv, equal_var=False)
        return t, p, len(av), len(bv), np.nan
    if mode == "paired_published":
        # Verbatim behaviour of the notebook cell: restrict both groups to the
        # culture numbers they share, then hand the two series to ttest_rel,
        # which pairs by row position and ignores the index that was just set.
        if keys_a is None or keys_b is None:
            av, bv = a, b
        else:
            common = set(a.index) & set(b.index)
            av = a[[i in common for i in a.index]]
            bv = b[[i in common for i in b.index]]
        if len(av) != len(bv) or len(av) < 2:
            return np.nan, np.nan, len(av), len(bv), np.nan
        t, p = stats.ttest_rel(av.values, bv.values)
        return t, p, len(av), len(bv), len(av)
    if mode == "paired_fixed":
        common = sorted(set(a.index) & set(b.index), key=lambda x: (len(str(x)), str(x)))
        av = a.reindex(common).astype(float)
        bv = b.reindex(common).astype(float)
        ok = av.notna() & bv.notna()
        av, bv = av[ok], bv[ok]
        if len(av) < 2:
            return np.nan, np.nan, len(a), len(b), len(av)
        t, p = stats.ttest_rel(av, bv)
        return t, p, len(a), len(b), len(av)
    raise ValueError(mode)


def panel_tests(gf, comparisons, value="IC50", group_col="group",
                modes=("welch", "paired_fixed", "paired_published")):
    """Reproduce historical August comparisons; not a current-analysis entry point.

    ``comparisons`` is a list of (group1, group2) tuples in the order the panel
    reports them.
    """
    gf = split_strain(gf) if "culturenumber" not in gf else gf
    rows = []
    for drug in sorted(gf["Antibiotic"].dropna().unique()):
        d = gf[gf["Antibiotic"] == drug]
        for g1, g2 in comparisons:
            s1 = d[d[group_col] == g1]
            s2 = d[d[group_col] == g2]
            row = dict(Antibiotic=drug, Group1=g1, Group2=g2, value=value)
            m1, h1, n1 = ci95(s1[value])
            m2, h2, n2 = ci95(s2[value])
            row.update(mean1=m1, ci95_1=h1, n1=n1, mean2=m2, ci95_2=h2, n2=n2)
            for mode in modes:
                t, p, _, _, npair = compare(
                    s1[value], s2[value], mode,
                    keys_a=s1["culturenumber"], keys_b=s2["culturenumber"])
                row[f"t_{mode}"] = t
                row[f"p_{mode}"] = p
                if mode == "paired_fixed":
                    row["n_pairs"] = npair
            # Superseded August selection, retained for historical reproduction.
            dm = design_mode(g1, g2)
            dm = dm if f"p_{dm}" in row else "welch"
            row["design_test"] = "paired" if dm == "paired_fixed" else "unpaired Welch"
            row["t_design"] = row.get(f"t_{dm}")
            row["p_design"] = row.get(f"p_{dm}")
            rows.append(row)
    return pd.DataFrame(rows)


def correct(df, pcol, family_cols=("Antibiotic",), methods=("bonferroni", "fdr_bh")):
    """Add corrected p-values within each family defined by ``family_cols``.

    Correction is applied over the comparisons a panel actually reports for one
    drug, with the stated comparison counts (k = 3 in Figure 2, k = 6 in
    Figure 3). NaN p-values are excluded from k, so a comparison that could not
    be computed does not inflate the penalty on the ones that could.
    """
    out = df.copy()
    for method in methods:
        out[f"{pcol}_{method}"] = np.nan
    out[f"{pcol}_k"] = np.nan
    for _, idx in out.groupby(list(family_cols)).groups.items():
        sub = out.loc[idx, pcol]
        ok = sub.notna()
        k = int(ok.sum())
        out.loc[idx, f"{pcol}_k"] = k
        if k == 0:
            continue
        for method in methods:
            adj = multipletests(sub[ok].values, alpha=0.05, method=method)[1]
            out.loc[sub[ok].index, f"{pcol}_{method}"] = adj
    return out


def stars(p):
    if pd.isna(p):
        return "n.d."
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def all_pairs(groups):
    return list(itertools.combinations(groups, 2))


# Preserve the superseded August all-Welch choice for historical callers only.
# This set does not describe the current experiment. See final_stats.py.
PAIRED_COHORTS = set()


def design_mode(group1, group2):
    """Historical August selector; current design decisions live in final_stats.py."""
    return ("paired_fixed"
            if group1 in PAIRED_COHORTS and group2 in PAIRED_COHORTS
            else "welch")
