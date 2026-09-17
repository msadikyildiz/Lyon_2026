"""Check the 36 revised primary results against the committed numerical reference."""
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'working/analysis/stats-rework'

def main():
    reference = pd.read_csv(SCRIPTS/'tests/primary_ic50_v2.csv')
    actual = pd.read_csv(SCRIPTS/'out/final_statistics.csv').query('value == "IC50"')
    keys = ['Panel', 'Antibiotic', 'Group1', 'Group2']
    actual = actual[reference.columns].sort_values(keys).reset_index(drop=True)
    reference = reference.sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(reference, actual, rtol=1e-8, atol=1e-10, check_exact=False)
    print('36 revised primary contrasts match the committed reference')

if __name__ == '__main__':
    main()
