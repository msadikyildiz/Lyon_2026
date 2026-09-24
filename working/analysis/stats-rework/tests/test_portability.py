"""Keep scientific failures distinct from platform-specific rendering changes."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from openpyxl import Workbook
from PIL import Image
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
from compare_outputs import compare, OUT, checked_single_cell_report
from compare_caches import compare_table
from plot_style import figure_font
import reproduce
import reproduce_genomic_plots as genomic_plots


class PortabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.reference = Path(self.temp.name) / 'reference'
        self.rebuilt = Path(self.temp.name) / 'rebuilt'
        for root in [self.reference, self.rebuilt]:
            for folder in [OUT/'figures-final', OUT/'print-panels', OUT/'trajectories',
                           OUT/'genomic-panels', Path('working/figures-assembled')]:
                (root/folder).mkdir(parents=True)
                Image.new('RGB', (12, 10), 'white').save(root/folder/'panel.png')
            (root/OUT/'values.csv').write_text('name,value\na,1.25\n', encoding='utf-8')
            for folder in ['data/genomics/generated_tables', 'data/genomics/plot_tables',
                           'data/single-cell/plot_tables']:
                (root/folder).mkdir(parents=True)
                (root/folder/'values.csv').write_text('name,value\na,1.25\n', encoding='utf-8')
            (root/'scripts').mkdir()
            script = root/'scripts/reproduce_genomic_plots.py'
            script.write_text('# fixture\n')
            source = 'data/genomics/generated_tables/values.csv'
            manifest = {'files': {source: hashlib.sha256((root/source).read_bytes()).hexdigest()},
                        'plot_script_sha256': hashlib.sha256(script.read_bytes()).hexdigest()}
            (root/'data/genomics/plot_tables/manifest.json').write_text(json.dumps(manifest))
            (root/'working/source-data').mkdir()
            book = Workbook(); book.active.append(['a', 1.25])
            book.save(root/'working/source-data/Source Data.xlsx')

    def test_font_has_bundled_fallback(self):
        with patch.dict('os.environ', {'LYON_PLOT_FONT': 'Absent Font 123456789'}):
            self.assertEqual(figure_font(), 'STIXGeneral')

    def test_fit_array_serialization_uses_numeric_tolerance(self):
        a = pd.DataFrame({'Strain': ['P1'], 'IC50': [1.0], 'ic50_bootstrap': ['[1.0, 2.0, nan]']})
        b = a.copy()
        b['ic50_bootstrap'] = ['[1.000000000001, 2, NaN]']
        compare_table(a, b)
        b['ic50_bootstrap'] = ['[1.01, 2, NaN]']
        with self.assertRaises(AssertionError):
            compare_table(a, b)

    def test_fit_array_missing_values_and_shapes_remain_strict(self):
        a = pd.DataFrame({'x_fit': ['[1.0, nan, 3.0]']})
        for value in ['[1.0, 2.0, 3.0]', '[1.0]', '[[1.0, nan, 3.0]]']:
            with self.assertRaises(AssertionError):
                compare_table(a, pd.DataFrame({'x_fit': [value]}))

    def test_fit_identity_remains_exact(self):
        a = pd.DataFrame({'Strain': ['P1'], 'x_fit': ['[1.0]']})
        b = a.copy()
        b['Strain'] = ['P2']
        with self.assertRaises(AssertionError):
            compare_table(a, b)

    def test_s9_order_is_independent_of_tied_input_rows(self):
        source = genomic_plots.table('combined_lineages')
        results = []
        for seed in [17, 83]:
            with patch.object(genomic_plots, 'OUT', Path(self.temp.name)), \
                 patch.object(genomic_plots, 'CROSSWALK', []), patch.object(genomic_plots, 'CHANGES', []), \
                 patch.object(genomic_plots, 'table', return_value=source.sample(frac=1, random_state=seed)):
                matrix, labels, _ = genomic_plots.s9()
                results.append((matrix, labels, list(genomic_plots.CROSSWALK)))
        pd.testing.assert_frame_equal(results[0][0], results[1][0], check_exact=True)
        self.assertEqual(results[0][1:], results[1][1:])

    def test_render_change_is_reported_and_strict_mode_rejects_it(self):
        Image.new('RGB', (12, 10), 'red').save(self.rebuilt/OUT/'figures-final/panel.png')
        report = compare(self.reference, self.rebuilt)
        self.assertEqual(report['status'], 'passed')
        self.assertEqual(len(report['pngs_different']), 1)
        with self.assertRaisesRegex(ValueError, 'PNG differs'):
            compare(self.reference, self.rebuilt, strict_images=True)

    def test_dimension_change_fails(self):
        Image.new('RGB', (11, 10)).save(self.rebuilt/'working/figures-assembled/panel.png')
        with self.assertRaisesRegex(ValueError, 'PNG dimensions differ'):
            compare(self.reference, self.rebuilt)

    def test_intermediate_crop_change_is_reported(self):
        Image.new('RGB', (11, 10)).save(self.rebuilt/OUT/'figures-final/panel.png')
        report = compare(self.reference, self.rebuilt)
        self.assertEqual(len(report['intermediate_crop_differences']), 1)
        with self.assertRaisesRegex(ValueError, 'PNG dimensions differ'):
            compare(self.reference, self.rebuilt, strict_images=True)

    def test_numeric_change_still_fails(self):
        (self.rebuilt/OUT/'values.csv').write_text('name,value\na,1.5\n')
        with self.assertRaises(AssertionError):
            compare(self.reference, self.rebuilt)

    def test_line_endings_do_not_change_results(self):
        (self.rebuilt/'data/single-cell/plot_tables/values.csv').write_bytes(b'name,value\r\na,1.25\r\n')
        self.assertEqual(compare(self.reference, self.rebuilt)['status'], 'passed')

    def test_single_cell_roundoff_uses_existing_numerical_tolerance(self):
        file = self.rebuilt/'data/single-cell/plot_tables/values.csv'
        file.write_text('name,value\na,1.2500000000000002\n')
        self.assertEqual(compare(self.reference, self.rebuilt)['status'], 'passed')
        file.write_text('name,value\na,1.2501\n')
        with self.assertRaises(AssertionError):
            compare(self.reference, self.rebuilt)

    def test_single_cell_validation_requires_small_errors_and_exact_counts(self):
        good = {'rows': 20, 'p_val_max_log10_difference': 1e-10}
        bad = {**good, 'p_val_max_log10_difference': 1e-5}
        with self.assertRaisesRegex(ValueError, 'exceeds tolerance'):
            checked_single_cell_report(bad)
        self.assertNotEqual(checked_single_cell_report(good),
                            checked_single_cell_report({**good, 'rows': 21}))

    def test_invalid_provenance_hash_fails(self):
        (self.rebuilt/'scripts/reproduce_genomic_plots.py').write_text('# changed\n')
        with self.assertRaises(AssertionError):
            compare(self.reference, self.rebuilt)

    def test_r_library_path_is_passed_without_shell_syntax(self):
        library = Path(self.temp.name)/'R library with spaces'
        with patch.object(reproduce, 'ROOT', self.rebuilt), patch.object(reproduce, 'run'), \
             patch.object(reproduce.platform, 'platform', return_value='test-platform'), \
             patch.object(reproduce.shutil, 'which', return_value='Rscript'), \
             patch.object(sys, 'argv', ['reproduce.py', '--single-cell', '--r-library', str(library)]), \
             patch.object(reproduce.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as call:
            reproduce.main()
        self.assertEqual(call.call_args.kwargs['env']['R_LIBS_USER'], str(library.resolve()))
        self.assertEqual(call.call_args.args[0][-1], '--dge')


if __name__ == '__main__':
    unittest.main()
