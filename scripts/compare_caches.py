"""Compare refitted values, dose-level records and all stored arrays with shipped caches."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ARRAY_COLUMNS = {'x_fit', 'y_fit', 'ic50_bootstrap'}
RTOL, ATOL = 1e-8, 1e-10


def compare_table(reference, actual):
    """Compare serialized fitted arrays numerically, retaining exact identities."""
    pd.testing.assert_index_equal(reference.index, actual.index)
    pd.testing.assert_index_equal(reference.columns, actual.columns)
    arrays = [col for col in reference if col in ARRAY_COLUMNS]
    pd.testing.assert_frame_equal(reference.drop(columns=arrays), actual.drop(columns=arrays),
                                  check_exact=False, rtol=RTOL, atol=ATOL)
    for column in arrays:
        for row, (a, b) in enumerate(zip(reference[column], actual[column])):
            a = json.loads(a.replace('nan', 'NaN')) if isinstance(a, str) else a
            b = json.loads(b.replace('nan', 'NaN')) if isinstance(b, str) else b
            a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
            if a.ndim != 1 or a.shape != b.shape:
                raise AssertionError(f'{column}, row {row}: fit array shapes differ or are not one-dimensional')
            np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL, equal_nan=True,
                                       err_msg=f'{column}, row {row}')


def main():
    reference, actual, output = map(Path, sys.argv[1:])
    expected = sorted(p.name for p in reference.glob('*.pkl'))
    if sorted(p.name for p in actual.glob('*.pkl')) != expected:
        raise ValueError('Cache file set differs')
    results = []
    for name in expected:
        a, b = pd.read_pickle(reference/name), pd.read_pickle(actual/name)
        compare_table(a, b)
        results.append({'file': name, 'rows': len(a), 'columns': len(a.columns), 'passed': True})
    output.write_text(json.dumps({'rtol': RTOL, 'atol': ATOL, 'tables': results}, indent=2)+'\n')
    print(f'{len(results)} cache tables match within rtol=1e-8, atol=1e-10')

if __name__ == '__main__':
    main()
