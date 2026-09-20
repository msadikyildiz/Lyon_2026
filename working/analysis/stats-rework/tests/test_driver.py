"""Regression tests for failed reruns and isolated refit output handling."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import reproduce as driver

class DriverTests(unittest.TestCase):
    def test_failure_replaces_previous_pass_record(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'runs').mkdir()
            report=root/'runs/reproduction.json';report.write_text('{"status":"passed"}')
            with patch.object(driver,'ROOT',root), patch.object(sys,'argv',['reproduce.py']), patch.object(driver,'run',side_effect=RuntimeError('fixture failure')):
                with self.assertRaisesRegex(RuntimeError,'fixture failure'):driver.main()
            self.assertEqual(json.loads(report.read_text())['status'],'failed')

    def test_refit_removes_only_copied_generated_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'source';source.mkdir()
            generated=[driver.SCRIPTS/'out',Path('working/source-data'),Path('working/figures-assembled'),
                       Path('data/genomics/generated_tables'),Path('data/genomics/plot_tables')]
            for folder in generated:
                (source/folder).mkdir(parents=True);(source/folder/'sentinel').write_text('reference')
            cache=source/driver.SCRIPTS/'cache';cache.mkdir();(cache/'fit.pkl').write_text('reference cache')
            output=Path(temp)/'refit'
            with patch.object(driver,'ROOT',source), patch.object(sys,'argv',['reproduce.py','--refit','--jobs','1','--output',str(output)]), patch.object(driver,'run',side_effect=RuntimeError('stop after setup')):
                with self.assertRaisesRegex(RuntimeError,'stop after setup'):driver.main()
            for folder in generated:
                self.assertTrue((source/folder/'sentinel').exists())
                self.assertFalse((output/folder).exists())
            self.assertEqual((output/'reference-cache/fit.pkl').read_text(),'reference cache')
            self.assertEqual(list((output/driver.SCRIPTS/'cache').iterdir()),[])

if __name__ == '__main__':unittest.main()
