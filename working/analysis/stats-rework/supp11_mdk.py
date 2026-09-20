"""Supplementary Figure 11 (double-mutant MDK) recomputed from colony counts.

Why: the workbook MDK_double_mutant.xlsx computes the time-zero CFU per mL
(count / dilution / 0.01 mL) but stores the 3 h and 7 h CFU as count / dilution
only. The confirmed plated volume was 10 uL at every measured time point.
Adam's Figure 5e workbook normalises every time point the same way and puts
hipA at 8.4e-5 (3 h) and 1.4e-5 (7 h); the double-mutant workbook has the same
strain at 1.3e-6 and 1.1e-7.

Fix: survivor fraction = (count_t / dilution_t) / (count_0 / dilution_0), with
the same 10 uL plated volume at every time point (confirmed 19 September 2026), so the volume
cancels. A count of zero is a detection-limit observation and is kept, plotted
at the fraction one colony would give, and flagged. Replicate spread is the
scaled MAD (scipy median_abs_deviation, scale='normal'), which is what both
MDK notebooks in the paper use.

Run:  python supp11_mdk.py   -> out/supp11_recomputed.csv, out/supp11_summary.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
XLSX = (HERE.parent.parent / "figures" / "Supplemental Figure 11 - double mutants"
        / "MDK_double_mutant.xlsx")
LABEL = {"wt": "MG", "hipA": "MG hipA", "hipA_selB": "MG hipA+selB",
         "hipA_ftsH": "MG hipA+ftsH", "hipA_rpoZ": "MG hipA+rpoZ"}


def main():
    d = pd.read_excel(XLSX, sheet_name="Sheet1")
    d['source_excel_row'] = d.index + 2
    d = d[d["Name"].notna()].copy()
    d["Replicate"] = d["Replicate"].astype(int)
    d["Time"] = d["Time"].astype(int)
    if d.duplicated(['Name', 'Replicate', 'Time']).any() or len(d) != 60:
        raise ValueError('Expected 60 unique strain/replicate/time observations')
    if not np.isfinite(d[['Count', 'Factor']]).all().all() or (d.Count < 0).any() or (d.Factor <= 0).any():
        raise ValueError('Missing/invalid count or dilution factor')
    if set(d.Name) != set(LABEL):
        raise ValueError('Unexpected strains')
    for _, g in d.groupby('Name'):
        if set(zip(g.Replicate, g.Time)) != {(r, t) for r in range(1, 5) for t in [0, 3, 7]}:
            raise ValueError('Missing replicate/time observation')
    # colonies per plated volume, corrected for the plate dilution
    d["cfu_per_plated_volume"] = d["Count"] / d["Factor"]
    t0 = d[d["Time"] == 0].set_index(["Name", "Replicate"])["cfu_per_plated_volume"]
    if (t0 <= 0).any():
        raise ValueError('Time-zero counts must be positive')
    d["cfu_time_zero"] = [t0[(n, r)] for n, r in zip(d["Name"], d["Replicate"])]
    d["detection_limit"] = (1 / d["Factor"]) / d["cfu_time_zero"]
    d["below_detection"] = d["Count"] == 0
    measured = d['cfu_per_plated_volume'] / d['cfu_time_zero']
    d['fraction'] = measured.where(~d.below_detection)
    d['fraction_lower'] = measured
    d['fraction_upper'] = np.where(d.below_detection, d.detection_limit, measured)
    d['fraction_plot'] = d.fraction_upper
    d['normalization_assumption'] = '10 uL at every measured time point; confirmed 19 September 2026'
    d['concentration_assumption'] = 'Factor follows recorded dilution, including concentration when Factor > 1'
    d['censoring'] = np.where(d.below_detection, 'below one-colony detection limit', 'measured')
    d["fraction_workbook"] = d["Fraction"]
    d["workbook_to_recomputed"] = d["fraction_workbook"] / d["fraction"]
    d["strain"] = d["Name"].map(LABEL)
    cols = ['source_excel_row', "strain", "Name", "Replicate", "Time", "Count", "Dilution", "Factor",
            "cfu_per_plated_volume", "cfu_time_zero", "fraction", 'fraction_lower', 'fraction_upper',
            'fraction_plot', "below_detection", 'censoring', "detection_limit", "fraction_workbook",
            "workbook_to_recomputed", 'normalization_assumption', 'concentration_assumption']

    s = (d.groupby(["Name", "Time"])["fraction_plot"]
         .agg(n="size", median="median",
              mad_scaled=lambda x: median_abs_deviation(x, scale="normal"),
              min="min", max="max")
         .reset_index())
    hip = s[s.Name == "hipA"].set_index("Time")["median"]
    s["ratio_to_hipA_median"] = [m / hip[t] if t in hip.index and t > 0 else np.nan
                                 for m, t in zip(s["median"], s["Time"])]
    # replicate-level separation from hipA: how many of the strain's replicates
    # exceed the highest hipA replicate at that time point
    hmax = d[d.Name == "hipA"].groupby("Time")["fraction"].max()
    s["replicates_above_all_hipA"] = [
        int((d[(d.Name == n) & (d.Time == t)]["fraction"] > hmax[t]).sum()) if t > 0 else np.nan
        for n, t in zip(s["Name"], s["Time"])]
    s["strain"] = s["Name"].map(LABEL)
    s["n_below_detection"] = [int(d[(d.Name == n) & (d.Time == t)]["below_detection"].sum())
                              for n, t in zip(s["Name"], s["Time"])]
    s['mad_scaled_min'] = s.mad_scaled
    s['mad_scaled_max'] = s.mad_scaled
    for idx, row in s[s.n_below_detection > 0].iterrows():
        g = d[(d.Name == row.Name) & (d.Time == row.Time)]
        # For this dataset the single censored value is below every observed value.
        # The median is fixed, and MAD is monotone over its allowed interval.
        if len(g) != 4 or row.n_below_detection != 1 or g.loc[g.below_detection, 'fraction_upper'].max() >= g.fraction.min():
            raise ValueError('Censoring configuration requires a new bound calculation')
        endpoints = [g.fraction_lower.values, g.fraction_upper.values]
        if not np.isclose(np.median(endpoints[0]), np.median(endpoints[1]), rtol=1e-12, atol=0):
            raise ValueError('Median depends on censoring substitution')
        bounds = [median_abs_deviation(v, scale='normal') for v in endpoints]
        s.loc[idx, ['mad_scaled_min', 'mad_scaled_max']] = [min(bounds), max(bounds)]
    s['spread_convention'] = 'normal-scaled MAD; central displayed value uses censored upper bound'
    s['normalization_assumption'] = d.normalization_assumption.iloc[0]
    OUT.mkdir(exist_ok=True)
    d[cols].to_csv(OUT / 'supp11_recomputed.csv', index=False)
    s.to_csv(OUT / "supp11_summary.csv", index=False)

    print("workbook / recomputed fraction at t > 0 (should be a constant if only the volume was dropped):")
    print(d[d.Time > 0]["workbook_to_recomputed"].describe()[["min", "50%", "max"]].to_string())
    print("\nmedian survivor fraction, scaled MAD, ratio to hipA, replicates above every hipA replicate:")
    print(s[s.Time > 0][["strain", "Time", "n", "median", "mad_scaled", "ratio_to_hipA_median",
                         "replicates_above_all_hipA", "n_below_detection"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
