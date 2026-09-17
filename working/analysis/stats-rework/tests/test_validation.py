"""Failure cases for lineage matching and publication-output validation."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import final_stats as f


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.a = pd.Series([1., 2., 4.], index=['1', '2', '3'])
        self.b = pd.Series([2., 3., 9.], index=['1', '2', '3'])

    def test_order_does_not_change_pairs(self):
        ordered = f.paired(self.a, self.b)
        shuffled = f.paired(self.a.iloc[::-1], self.b.iloc[[1, 2, 0]])
        for key in ['ratio', 'lo', 'hi', 'p', 't']:
            self.assertEqual(ordered[key], shuffled[key])
        self.assertEqual(shuffled['_keys'], ['1', '2', '3'])

    def test_missing_pair_is_not_silently_dropped(self):
        with self.assertRaisesRegex(ValueError, 'Incomplete pairs'):
            f.paired(self.a, self.b.iloc[:2])

    def test_duplicate_ids_fail(self):
        b = self.b.copy()
        b.index = ['1', '2', '2']
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            f.paired(self.a, b)

    def test_invalid_measurements_fail(self):
        for value in [0., -1., np.nan, np.inf]:
            b = self.b.copy()
            b.iloc[0] = value
            with self.assertRaisesRegex(ValueError, 'finite positive'):
                f.paired(self.a, b)

    def test_empty_paired_draw_array_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Empty or invalid'):
            f.nested_boot([np.array([0., 1.]), np.array([])],
                          [np.array([0., 1.]), np.array([1., 2.])], 'paired')

    def test_failed_validation_preserves_existing_outputs(self):
        # An exclusion failure occurs before any numerical calculation or CSV replacement.
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            sentinel = out / 'final_statistics.csv'
            sentinel.write_bytes(b'previous valid results\n')
            with patch.object(f, 'OUT', out), patch.object(f, 'load', return_value=pd.DataFrame()), \
                 patch.object(f, 'validate_inputs', return_value=[]), \
                 patch.object(f, 'exclusion_checks', return_value=[('bad exclusion', '', '', '', False)]), \
                 patch.object(f, 'contrast_rows') as contrast:
                with self.assertRaisesRegex(ValueError, 'Input validation failed'):
                    f.main()
                contrast.assert_not_called()
            self.assertEqual(sentinel.read_bytes(), b'previous valid results\n')
            self.assertEqual(list(out.iterdir()), [sentinel])


if __name__ == '__main__':
    unittest.main()
