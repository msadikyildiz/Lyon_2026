"""Regenerate growth features for every panel that carries a statistical test.

Writes one pickle per dataset into ``cache/`` so the curve fits are paid for
once. Each entry mirrors one notebook: same input workbook, same Group filter,
same thresholds.

Inputs come from the BioHPC copy under ``biohpc-pull/data``. Those are the exact
files that produced the published numbers. The OneDrive copies carry identical
values but ``ODFinal.xlsx`` there has already had its stray trailing column
removed, so the notebooks' verbatim ``.iloc[:, :-1]`` would delete a real dose
column if pointed at them. Verified: local ODFinal equals BioHPC ODFinal's first
14 columns, and both PlateInfo workbooks are value-identical to their local
twins.

Run:  python build_all.py [dataset ...]
"""

import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from pipeline import (ATEC_TO_REMOVE, ATEC_TO_REMOVE_ALL_DRUGS,
                      GROUPS_TO_REMOVE_FULL, GROUPS_TO_REMOVE_SHORT,
                      PLATE_ID_OFFSET, PLATE_ID_OFFSET_UNT,
                      STRAINS_TO_REMOVE_FULL, STRAINS_TO_REMOVE_SHORT,
                      build_growth_features)

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
CACHE.mkdir(exist_ok=True)
DATA = HERE / "biohpc-pull" / "data"

PAPLPC = (DATA / "ODFinal.xlsx", DATA / "PlateInfo.xlsx")
MUTANTS = (DATA / "ODFinal_mutants.xlsx", DATA / "PlateInfo_mutants_norep.xlsx")
UNT = (DATA / "ODFinal_unt.xlsx", DATA / "PlateInfo_unt.xlsx")

FULL = dict(groups_to_remove=GROUPS_TO_REMOVE_FULL,
            strains_to_remove=STRAINS_TO_REMOVE_FULL,
            drug_strains_to_remove=None)
# Short lists plus the cefepime-only ATEC filter: what Figure 3, Supplementary
# Figure 3, Supplementary Figure 4 and Supplementary Figure 7 actually ran.
SHORT = dict(groups_to_remove=GROUPS_TO_REMOVE_SHORT,
             strains_to_remove=STRAINS_TO_REMOVE_SHORT,
             drug_strains_to_remove=ATEC_TO_REMOVE)
NONE = dict(groups_to_remove=[], strains_to_remove=[],
            drug_strains_to_remove=None)

DATASETS = {
    "fig2_paplpc":  (*PAPLPC,  ["P", "PA", "PC", "PL"],
                     "Figure 2b-d, Supplementary Figure 2",
                     dict(plate_offsets=PLATE_ID_OFFSET, **FULL)),
    "fig3_plac":    (*PAPLPC,  ["P", "PL", "PLA", "PLAC"],
                     "Figure 3b-d, Supplementary Figure 7",
                     dict(plate_offsets=PLATE_ID_OFFSET, **SHORT)),
    "supp3_pcr":    (*PAPLPC,  ["P", "PCr"],
                     "Supplementary Figure 3b-d",
                     dict(plate_offsets=PLATE_ID_OFFSET, **SHORT)),
    # As published: plating-error cultures dropped for cefepime only.
    "supp4_atec":   (*PAPLPC,  ["ATEC", "ATEC-C", "ATEC-C-r"],
                     "Supplementary Figure 4b-d (as published)",
                     dict(plate_offsets=PLATE_ID_OFFSET, **SHORT)),
    # Counterfactual: the same exclusions applied to all three drugs, so the
    # effect of the drug-specific filter can be quantified rather than asserted.
    "supp4_atec_excluded": (*PAPLPC, ["ATEC", "ATEC-C", "ATEC-C-r"],
                     "Supplementary Figure 4b-d (plating errors excluded for all drugs)",
                     dict(plate_offsets=PLATE_ID_OFFSET,
                          groups_to_remove=GROUPS_TO_REMOVE_SHORT,
                          strains_to_remove=STRAINS_TO_REMOVE_SHORT,
                          drug_strains_to_remove=ATEC_TO_REMOVE_ALL_DRUGS)),
    "fig5_mutants": (*MUTANTS, None,
                     "Figure 5b-d, Supplementary Figure 10a-c",
                     dict(plate_offsets=PLATE_ID_OFFSET, **NONE)),
    "supp6_unt":    (*UNT,     ["P", "P-unt"],
                     "Supplementary Figure 6a-c",
                     dict(plate_offsets=PLATE_ID_OFFSET_UNT, **NONE)),
}


DATASETS.pop("supp4_atec_excluded")  # Withdrawn counterfactual is not a supported entry point.

def main(names):
    for name in names:
        od, pi, groups, panels, kw = DATASETS[name]
        out = CACHE / f"{name}.pkl"
        if out.exists():
            print(f"[skip] {name} already cached", flush=True)
            continue
        t0 = time.time()
        print(f"[run ] {name}  ({panels})", flush=True)
        gf, dfa = build_growth_features(od, pi, groups, **kw)
        gf.to_pickle(out)
        dfa.to_pickle(CACHE / f"{name}__df_analysis.pkl")
        print(f"[done] {name}  growth_features={gf.shape}  {time.time() - t0:.0f}s",
              flush=True)
    print("[ALL DONE]", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or [name for name in DATASETS if name != "supp4_atec_excluded"])
