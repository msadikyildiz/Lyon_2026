"""Dose-response fitting with the pinned plategig pipeline.

Reproduces cells 1-26 of ``241011_Adam_mic_*.ipynb`` so growth features can be
regenerated outside Jupyter. Historical numerical agreement is checked by ``validate.py`` against the saved notebook values.

plategig is vendored at commit 88839a2 to preserve the original fitting API.
"""

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.home() / ".cache/matplotlib"))
# One BLAS thread per process: the fits are tiny, and oversubscription costs
# more than it buys once several datasets run side by side.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import matplotlib
matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "plategig-88839a2"))

import pandas as pd
import plategig  # noqa: E402  (must follow the sys.path insert)

PLATEGIG_COMMIT = "88839a2"


class _SerialPool:
    """Drop-in for ``multiprocessing.Pool`` that runs starmap in this process.

    ``robust_phenotyper`` opens a fresh 8-process pool for every strain-drug
    combination, 1000 bootstrap fits each. Python on macOS spawns rather than
    forks, so each pool starts eight interpreters that re-import dependencies.
    Serial inner execution avoids this process-startup overhead.

    This changes no result. ``bootstrap_resample`` is called as
    ``[(x, y, p0, seed) for seed in range(n_bootstrap)]``, so every draw is
    seeded by its index and the outputs are identical in any execution order.
    Parallelism is applied one level up instead, across datasets.
    """

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def starmap(self, fn, iterable):
        return [fn(*args) for args in iterable]

    def close(self):
        pass

    def join(self):
        pass


plategig.plategig.mp.Pool = _SerialPool

# Experiment -> Plate_ID offset, from each notebook's calculate_plate_id. The
# mapping is per-notebook: the untreated workbook holds only Experiment 4 and
# maps it to no offset, where the PAPLPC workbook maps Experiment 1 that way.
PLATE_ID_OFFSET = {1: 0, 2: 99, 3: 148}
PLATE_ID_OFFSET_UNT = {4: 0, 2: 99, 3: 148}

# Cell 18 exclusions, keyed on the string Experiment column.
# Recorded reasons:
#   PAC exp 1        not used
#   PCr exp 1        dose range too small
#   PLAC7, PLAC10    bad curve in exp 1, remeasured in exp 2
#   ATEC 1, 2, 4     plating error
#   ATEC-C3          plating error
#   ATEC-C-R         plating error
#
# These lists are NOT the same in every notebook, and the difference matters in
# exactly one place.
#
#   Figure 2, Supplementary Figure 2   full group and strain lists, no drug-
#                                      specific filter. The ATEC entries are
#                                      inert because ATEC is not in the panel.
#   Figure 3, Supplementary 3, 4, 7    short group and strain lists plus
#                                      ATEC_TO_REMOVE, which drops the
#                                      plating-error cultures for CEFEPIME ONLY.
#   Supplementary 5a, 6a, Figure 5     no exclusions at all.
#
# Experiment-2 plates 39-41 are cefepime-only, so exclusions are drug-specific.
# The experimental account attributes the failed ATEC-C3 technical series to
# inoculation error; the available screenshot alone does not establish cause.
#
# Note also that ('2', 'ATEC-C-R', 'Cefepime') never matches anything: the group
# is spelled 'ATEC-C-r' in the data, and the filter is applied to Strain, which
# by this point carries a culture number appended.
GROUPS_TO_REMOVE_FULL = [("1", "PAC"), ("1", "PCr"), ("2", "ATEC-C-R")]
STRAINS_TO_REMOVE_FULL = [
    ("1", "PLAC7"), ("1", "PLAC10"),
    ("2", "ATEC1"), ("2", "ATEC2"), ("2", "ATEC4"), ("2", "ATEC-C3"),
]
GROUPS_TO_REMOVE_SHORT = [("1", "PAC"), ("1", "PCr")]
STRAINS_TO_REMOVE_SHORT = [("1", "PLAC7"), ("1", "PLAC10")]

# (Experiment, Strain, Antibiotic) triples, applied after the two filters above.
# 15 September 2026: the first entry was ("2", "ATEC-C-R", "Cefepime") in the
# notebook and matched nothing, because the group is spelled ATEC-C-r and the
# filter runs on Strain, which carries the culture number. The intended key is
# ATEC-C-r1. With the inert key, experiment-2 and experiment-3 cefepime curves
# for ATEC-C-r were pooled into one fit; the corrected key drops experiment 2
# as intended. ATEC-C-r is a single culture and enters no statistical contrast.
ATEC_TO_REMOVE = [
    ("2", "ATEC-C-r1", "Cefepime"), ("2", "ATEC1", "Cefepime"),
    ("2", "ATEC2", "Cefepime"), ("2", "ATEC4", "Cefepime"),
    ("2", "ATEC-C3", "Cefepime"),
]
def _plate_id(row, offsets=PLATE_ID_OFFSET):
    exp = row["Experiment"]
    return row["Plate"] + offsets[exp] if exp in offsets else None


def build_growth_features(od_path, plate_info_path, groups,
                          ic50_threshold=0.5, mic_threshold=0.05,
                          plate_offsets=None,
                          groups_to_remove=None, strains_to_remove=None,
                          drug_strains_to_remove=None):
    """Run the published pipeline and return (growth_features, df_analysis).

    ``groups`` is the Group whitelist applied at cell 19; pass None to keep all.
    ``plate_offsets``, ``groups_to_remove`` and ``strains_to_remove`` default to
    the PAPLPC notebook's values and should be set to whatever the notebook
    behind the panel actually used.
    """
    plate_offsets = PLATE_ID_OFFSET if plate_offsets is None else plate_offsets
    groups_to_remove = (GROUPS_TO_REMOVE_FULL if groups_to_remove is None
                        else groups_to_remove)
    strains_to_remove = (STRAINS_TO_REMOVE_FULL if strains_to_remove is None
                         else strains_to_remove)

    df_plate_info = pd.read_excel(plate_info_path, engine="openpyxl")
    df_plate_info["Plate_ID"] = df_plate_info.apply(
        lambda r: _plate_id(r, plate_offsets), axis=1)

    df_od_raw = pd.read_excel(od_path)
    df_od_raw = df_od_raw.iloc[:, :-1]  # trailing all-NaN column, an input mistake

    df = plategig.static.convert_OD_plate_to_long(df_od_raw, df_plate_info)

    background = plategig.static.calc_median_background_all_plates(
        df, df_plate_info, plot=False)

    df_bc = df.copy()
    df_bc["OD_final"] = df_bc["OD"] - background

    df_analysis = pd.merge(df_plate_info, df_bc, on=["Plate_ID", "Well"])
    df_analysis = df_analysis[[
        "Experiment", "Strain", "Culture", "Replicate", "Antibiotic",
        "Dose", "Plate_ID", "Well", "Row", "Column", "OD", "OD_final",
    ]]
    df_analysis = df_analysis[~df_analysis["Strain"].isin(["Media Only", "Cells Only"])]
    df_analysis["Group"] = df_analysis.Strain
    df_analysis["Strain"] = (
        df_analysis.Strain + df_analysis.Culture.astype(int).astype(str))

    df_analysis["Experiment"] = df_analysis["Experiment"].astype(str).str.strip()
    df_analysis = df_analysis[
        ~df_analysis[["Experiment", "Group"]].apply(tuple, axis=1).isin(groups_to_remove)]
    df_analysis = df_analysis[
        ~df_analysis[["Experiment", "Strain"]].apply(tuple, axis=1).isin(strains_to_remove)]
    if drug_strains_to_remove:
        df_analysis = df_analysis[
            ~df_analysis[["Experiment", "Strain", "Antibiotic"]]
            .apply(tuple, axis=1).isin(drug_strains_to_remove)]

    if groups is not None:
        df_analysis = df_analysis[df_analysis["Group"].isin(groups)]

    valid_combinations = plategig.static.prep_valid_combinations(
        df_analysis,
        multiplex=["Strain", "Antibiotic"],
        ic50_threshold=ic50_threshold,
        mic_threshold=mic_threshold)

    growth_features = plategig.static.apply_phenotyper(df_analysis, valid_combinations)
    growth_features = plategig.static.cap_growth_features_within_experiment_range(
        growth_features)

    return growth_features, df_analysis
