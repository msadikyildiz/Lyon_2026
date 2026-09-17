"""Compare refitted values, dose-level records and all stored arrays with shipped caches."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd


def main():
    reference, actual, output = map(Path, sys.argv[1:])
    expected = sorted(p.name for p in reference.glob('*.pkl'))
    if sorted(p.name for p in actual.glob('*.pkl')) != expected:
        raise ValueError('Cache file set differs')
    results = []
    for name in expected:
        a, b = pd.read_pickle(reference/name), pd.read_pickle(actual/name)
        # pandas also compares arrays nested in object columns. NaNs must occupy the same cells.
        pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=1e-8, atol=1e-10)
        results.append({'file': name, 'rows': len(a), 'columns': len(a.columns), 'passed': True})
    output.write_text(json.dumps({'rtol': 1e-8, 'atol': 1e-10, 'tables': results}, indent=2)+'\n')
    print(f'{len(results)} cache tables match within rtol=1e-8, atol=1e-10')

if __name__ == '__main__':
    main()
