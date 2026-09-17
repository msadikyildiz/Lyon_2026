"""Growth-data loader for the descriptive growth analysis.

Adam confirmed six technical wells from one culture per strain (biological
n=1). Use final_stats.py for current summaries. The inferential main routine
below is historical and is disabled; its earlier output remains in Git.
Doubling times are minutes despite the original notebook's hour label.
"""

import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import linregress


HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
XLSX = (HERE.parent.parent / "figures" /
        "Figure 5 - mutants" / "A - doubling times" / "growth_022525_mutants.xlsx")

WELL_CONTENTS = {
    "gatA": [f"{r}3" for r in "BCDEFG"],
    "glvC": [f"{r}4" for r in "BCDEFG"],
    "selB": [f"{r}5" for r in "BCDEFG"],
    "hipA": [f"{r}6" for r in "BCDEFG"],
    "rpoZ": [f"{r}7" for r in "BCDEFG"],
    "ftsH": [f"{r}8" for r in "BCDEFG"],
    "fimE": [f"{r}9" for r in "BCDEFG"],
    "MG":   [f"{r}10" for r in "BCDEFG"],
}
OD_MIN, OD_MAX = 0.02, 0.08


def simple_growth_parameters(time, od, od_min=OD_MIN, od_max=OD_MAX):
    """Verbatim from growth.ipynb cell 2."""
    valid = (od > 0) & (od >= od_min) & (od <= od_max)
    t, o = time[valid], od[valid]
    if len(t) < 2:
        return None, None, None, None
    slope, intercept, _, _, _ = linregress(t, np.log10(o))
    growth_rate = slope * np.log(10)
    doubling_time = (np.log(2) / growth_rate) * 60  # minutes, despite the "(h)" label
    return growth_rate, doubling_time, slope, intercept


def load_doubling_times():
    data = pd.read_excel(XLSX, sheet_name="Sheet1", header=30)
    data.drop(columns=["Unnamed: 0", "T° od600:600"], inplace=True)
    data["Time"] = pd.to_datetime(data["Time"], format="%H:%M:%S", errors="coerce")
    data["Time"] = (data["Time"].dt.hour + data["Time"].dt.minute / 60
                    + data["Time"].dt.second / 3600)
    data.dropna(subset=["Time"], inplace=True)
    data = data.loc[:data.dropna(how="all",
                                 subset=data.columns.difference(["Time"])).index[-1]]

    media_wells = ([f"A{c}" for c in range(1, 13)] + [f"H{c}" for c in range(1, 13)]
                   + [f"{r}{c}" for r in "BCDEFG" for c in [1, 2, 11, 12]])
    media_control = data[media_wells].mean(axis=1)
    data_bc = data.copy()
    data_bc.iloc[:, 1:] = (data.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
                           .subtract(media_control, axis=0))

    filt = data_bc[data_bc["Time"] <= 6]
    rows = []
    for culture, wells in WELL_CONTENTS.items():
        for well in wells:
            od = pd.to_numeric(filt[well], errors="coerce").values
            keep = ~np.isnan(od)
            t, o = filt["Time"].values[keep][2:], od[keep][2:]
            gr, dt, _, _ = simple_growth_parameters(t, o)
            rows.append(dict(Culture=culture, Well=well,
                             growth_rate_per_h=gr, doubling_time_min=dt))
    return pd.DataFrame(rows)


