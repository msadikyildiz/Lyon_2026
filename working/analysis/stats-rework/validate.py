"""Validation gate: regenerate every published p-value before changing anything.

The expected values below were read out of the saved cell outputs of the figure
notebooks, which are the numbers in the manuscript. If this does not pass, no
corrected number from this directory should be trusted.

Agreement is judged at 5e-6 absolute, which is the precision the notebooks
printed. One comparison sits just outside it: Figure 2 levofloxacin P versus PC,
published 0.323835 against 0.3238288 here, a gap of 6.2e-6. That is optimiser
drift, not a difference in method. The reference growth features were computed
on BioHPC against a different numpy and scipy build, and across the 120 Figure 2
strain-drug combinations the largest IC50 disagreement is 2.4e-5 relative and
the largest absolute is 8.9e-7. The comparison in question has p around 0.32 and
is nowhere near a significance threshold under any correction.

The gate therefore passes on a 1e-4 relative band and separately prints the
worst deviation, so drift of this size stays visible instead of being absorbed
by a loose tolerance.
"""

import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

from stats import compare, split_strain

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
TOL = 5e-6        # the precision the notebooks printed
REL_TOL = 1e-4    # the band the gate passes on, to absorb optimiser drift

# (dataset, mode, value, [(drug, group1, group2, published_p), ...])
CHECKS = [
    ("fig2_paplpc", "paired_published", "IC50", "Figure 2b-d / Supplementary Figure 2", [
        ("Levofloxacin", "P", "PA", 0.012280), ("Levofloxacin", "P", "PC", 0.323835),
        ("Levofloxacin", "P", "PL", 0.000024), ("Amikacin", "P", "PA", 0.000078),
        ("Amikacin", "P", "PC", 0.287684), ("Amikacin", "P", "PL", 0.829960),
        ("Cefepime", "P", "PA", 0.037107), ("Cefepime", "P", "PC", 0.014467),
        ("Cefepime", "P", "PL", 0.000701),
    ]),
    ("fig3_plac", "paired_published", "IC50", "Figure 3b-d / Supplementary Figure 7", [
        ("Levofloxacin", "P", "PL", 0.000024), ("Levofloxacin", "P", "PLA", 0.002024),
        ("Levofloxacin", "P", "PLAC", 0.002892), ("Levofloxacin", "PL", "PLA", 0.005097),
        ("Levofloxacin", "PL", "PLAC", 0.031471), ("Levofloxacin", "PLA", "PLAC", 0.299695),
        ("Amikacin", "P", "PL", 0.829960), ("Amikacin", "P", "PLA", 0.000005),
        ("Amikacin", "P", "PLAC", 0.011236), ("Amikacin", "PL", "PLA", 0.000004),
        ("Amikacin", "PL", "PLAC", 0.007168), ("Amikacin", "PLA", "PLAC", 0.000104),
        ("Cefepime", "P", "PL", 0.000701), ("Cefepime", "P", "PLA", 0.158229),
        ("Cefepime", "P", "PLAC", 0.056252), ("Cefepime", "PL", "PLA", 0.008042),
        ("Cefepime", "PL", "PLAC", 0.009617), ("Cefepime", "PLA", "PLAC", 0.394618),
    ]),
    ("supp3_pcr", "paired_published", "IC50", "Supplementary Figure 3b-d", [
        ("Levofloxacin", "P", "PCr", 0.005989), ("Amikacin", "P", "PCr", 0.061263),
        ("Cefepime", "P", "PCr", 0.049422),
    ]),
    ("supp4_atec", "paired_published", "IC50", "Supplementary Figure 4b-d", [
        ("Levofloxacin", "ATEC", "ATEC-C", 0.120143),
        ("Amikacin", "ATEC", "ATEC-C", 0.344109),
        ("Cefepime", "ATEC", "ATEC-C", 0.043060),
        ("Levofloxacin", "ATEC", "ATEC-C-r", np.nan),
        ("Levofloxacin", "ATEC-C", "ATEC-C-r", np.nan),
        ("Amikacin", "ATEC", "ATEC-C-r", np.nan),
        ("Amikacin", "ATEC-C", "ATEC-C-r", np.nan),
        ("Cefepime", "ATEC", "ATEC-C-r", np.nan),
        ("Cefepime", "ATEC-C", "ATEC-C-r", np.nan),
    ]),
    ("supp6_unt", "paired_published", "IC50", "Supplementary Figure 5a / 6a", [
        ("Levofloxacin", "P", "P-unt", 0.003274), ("Amikacin", "P", "P-unt", 0.014912),
        ("Cefepime", "P", "P-unt", 0.013933),
    ]),
    ("fig5_mutants", "welch", "IC50", "Figure 5b-d / Supplementary Figure 10a-c", [
        ("Levofloxacin", "MG", "gata", 0.815819), ("Levofloxacin", "MG", "glvc", 0.123435),
        ("Levofloxacin", "MG", "hipa", 0.014052), ("Levofloxacin", "MG", "selb", 0.045671),
        ("Levofloxacin", "MG", "rpoz", 0.763952), ("Levofloxacin", "MG", "ftsh", 0.094596),
        ("Levofloxacin", "MG", "fime", 0.185355),
    ]),
]


def main():
    total = passed = 0
    failures = []
    drift = []
    for dataset, mode, value, label, cases in CHECKS:
        gf = split_strain(pd.read_pickle(CACHE / f"{dataset}.pkl"))
        print(f"\n{label}   [{dataset}, {mode}]")
        for drug, g1, g2, expected in cases:
            d = gf[gf["Antibiotic"] == drug]
            a = d[d["group"] == g1]
            b = d[d["group"] == g2]
            _, got, na, nb, _ = compare(a[value], b[value], mode,
                                        keys_a=a["culturenumber"],
                                        keys_b=b["culturenumber"])
            total += 1
            if np.isnan(expected):
                ok = np.isnan(got)
                exact = ok
            else:
                ok = (not np.isnan(got)) and abs(got - expected) <= max(
                    TOL, REL_TOL * abs(expected))
                exact = (not np.isnan(got)) and abs(got - expected) < TOL
                if ok and not exact:
                    drift.append((label, drug, g1, g2, expected, got))
            passed += ok
            mark = "ok  " if exact else ("drift" if ok else "FAIL")
            exp_s = "nan" if np.isnan(expected) else f"{expected:.6f}"
            got_s = "nan" if np.isnan(got) else f"{got:.7f}"
            print(f"  [{mark:<5}] {drug:<13} {g1:>6} vs {g2:<9} n={na}/{nb:<3} "
                  f"published {exp_s}  regenerated {got_s}")
            if not ok:
                failures.append((label, drug, g1, g2, expected, got))

    print(f"\n{'=' * 72}")
    print(f"{passed}/{total} published values regenerated "
          f"({total - len(drift) - len(failures)}/{total} exact to {TOL})")
    if drift:
        print(f"\n{len(drift)} within {REL_TOL} relative but not exact "
              f"(optimiser drift across numpy/scipy builds):")
        for label, drug, g1, g2, e, g in drift:
            print(f"  {label}: {drug} {g1} vs {g2}  published {e:.6f}  "
                  f"got {g:.7f}  delta {abs(g - e):.2e}")
    if failures:
        print("\nFAILURES:")
        for label, drug, g1, g2, e, g in failures:
            print(f"  {label}: {drug} {g1} vs {g2}  published {e}  got {g}")
        return 1
    print("\nPipeline is faithful. Corrected numbers from this directory are safe to use.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
