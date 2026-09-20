"""Regenerate genomic tables from the supplied mutation calls using notebook cells."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path.home() / '.cache/matplotlib'))
os.environ['MPLBACKEND'] = 'Agg'
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/genomics'
DEST = DATA / 'generated_tables'
# Execute the original data-processing cells, including their annotation corrections.
# Display-only and plotting cells are excluded. Cell numbers are zero-based.
CONFIG = [
    ('MG_AMI', 'Figure 2/F-H heatmaps/05-PA.ipynb', [0, 1, 3, 5, 6, 9], 43),
    ('MG_LEV', 'Figure 2/F-H heatmaps/07-PL.ipynb', [0, 1, 3, 5, 6, 8, 9], 64),
    ('MG_CEF', 'Figure 2/F-H heatmaps/06-PC.ipynb', [0, 1, 3, 5, 6, 8, 12], 80),
    ('MG_CEF_R', 'Supplemental Figure 3 - CefR/E - genetics/08-PCr.ipynb', [0, 1, 3, 4, 5, 6, 7, 8], 56),
    ('PbEc_CEF', 'Supplemental Figure 4 - Pb/F - Mutations/10-ATEC12.ipynb', [0, 1, 4, 7, 8, 12, 14], 83),
    ('MG_untreated', 'Supplemental Figure 6/B - heatmap/05-Punt.ipynb', [0, 1, 4, 5, 6, 9], 38),
    ('PLAC', 'Figure 4/02-PLAC.ipynb', [0, 1, 4, 6, 9, 13, 14, 27], 601),
    ('combined_lineages', 'Supplemental Figure 9 - big heatmap/13-bigheatmap.ipynb', [0, 1, 3, 4, 5, 6, 8], 963),
]
ALIASES = {
    'Supplemental Figure 5 - Extended heatmaps/A/05-PA.ipynb': 'MG_AMI',
    'Supplemental Figure 5 - Extended heatmaps/B/07-PL.ipynb': 'MG_LEV',
    'Supplemental Figure 5 - Extended heatmaps/C/06-PC.ipynb': 'MG_CEF',
    'Supplemental Figure 8 - extended traced alleles/02-PLAC.ipynb': 'PLAC',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_inputs():
    manifest = json.loads((DATA / 'input_manifest.json').read_text())
    for item in manifest['files']:
        path = DATA / item['path']
        if path.stat().st_size != item['bytes'] or sha(path) != item['sha256']:
            raise ValueError(f'Input differs from supplied archive: {path}')
    meta = pd.read_csv(DATA / 'reference/metadata_complete.csv')
    required = {f"out/{r.FolderDate}/{r.source_file}/output/output.gd.tsv" for r in meta.itertuples()}
    supplied = {r['path'] for r in manifest['files'] if r['path'].endswith('.tsv')}
    if required != supplied:
        raise ValueError(f'Metadata/input mismatch: missing={required-supplied}, extra={supplied-required}')
    return manifest


def vector(value):
    if isinstance(value, str):
        # Notebook CSVs use both NumPy's space-separated and Python's comma-separated lists.
        return [float(x) for x in value.strip('[]').replace(',', ' ').split()]
    return list(value)


def scalar(value):
    if isinstance(value, (np.ndarray, list)):
        return list(value)
    if pd.isna(value):
        return None
    return value


def equivalent(a, b, displayed=False):
    a, b = scalar(a), scalar(b)
    if a is None or (isinstance(a, str) and a == ''):
        return b is None or (isinstance(b, str) and b == '')
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and bool(np.allclose(a, b, rtol=0, atol=1e-12))
    try:
        return abs(float(a) - float(b)) < 1e-12
    except (TypeError, ValueError):
        pass
    if displayed and isinstance(a, str) and a.endswith('...') and isinstance(b, str):
        return b.startswith(a[:-3])
    return a == b


def compare_table(reference, generated, keys, displayed=False):
    a, b = reference.copy(), generated.copy()
    for column in keys:
        for table in (a, b):
            if column == 'Pop':
                table[column] = table[column].astype(int)
            else:
                table[column] = table[column].fillna('').astype(str)
    if a.duplicated(keys).any() or b.duplicated(keys).any():
        raise ValueError(f'Nonunique mutation identity: {keys}')
    a, b = a.set_index(keys), b.set_index(keys)
    if not displayed and set(a.index) != set(b.index):
        raise ValueError('Mutation identities differ from supplied processed export')
    columns = [c for c in a.columns if c not in ['index', 'Unnamed: 0']]
    for key, row in a.iterrows():
        if key not in b.index:
            raise ValueError(f'Missing reference mutation: {key}')
        new = b.loc[key]
        for column in columns:
            av, bv = row[column], new[column]
            if column == 'freq':
                av, bv = vector(av), vector(bv)
            if not equivalent(av, bv, displayed):
                raise ValueError(f'Reference mismatch at {key}, {column}: {av!r} != {bv!r}')
    return len(a)


def regenerate(name, relative, cells, expected):
    path = ROOT / 'working/figures' / relative
    notebook = json.loads(path.read_text())
    old_cwd = Path.cwd()
    with tempfile.TemporaryDirectory(prefix='lyon-genomics-') as temporary:
        work = Path(temporary)
        (work / 'reproduce.py').touch()
        base = work / 'data/genomics'
        base.mkdir(parents=True)
        for folder in ['reference', 'out']:
            (base / folder).symlink_to(DATA / folder, target_is_directory=True)
        for lineage in ['PA', 'PL', 'PC', 'PCr', 'ATEC-C', 'PLAC']:
            (base / 'data/processed/traced_alleles' / lineage).mkdir(parents=True)
        namespace = {'__name__': '__main__'}
        try:
            os.chdir(work)
            for cell in cells:
                source = ''.join(notebook['cells'][cell]['source'])
                exec(compile(source, f'{relative}:cell-{cell}', 'exec'), namespace)
            table = namespace['combined_df'].copy()
            if len(table) != expected:
                raise ValueError(f'{name}: expected {expected} rows, obtained {len(table)}')
            # Preserve notebook calculations, while making exported row order deterministic.
            keys = [c for c in ['Pop', 'population', 'abs_position', 'position', 'ref_seq', 'new_seq', 'aa_new_seq'] if c in table]
            table = table.sort_values(keys, kind='stable').reset_index(drop=True)
            processed = []
            for generated in sorted((base / 'data/processed').rglob('*.csv')):
                supplied = DATA / generated.relative_to(base)
                if supplied.exists():
                    a, b = pd.read_csv(supplied), pd.read_csv(generated)
                    keys = ['abs_position'] if 'abs_position' in a else ['position', 'ref_seq', 'new_seq', 'aa_new_seq']
                    count = compare_table(a, b, keys)
                    processed.append({'file': str(supplied.relative_to(DATA)), 'rows_matched': count})
            return table, processed
        finally:
            os.chdir(old_cwd)
            plt.close('all')


def main():
    status = ROOT / 'runs/genomics.json'
    status.parent.mkdir(exist_ok=True)
    status.write_text('{"status": "running"}\n')
    try:
        manifest = verify_inputs()
        tables, entries, processed, comparisons = {}, {}, [], []
        recovered = json.loads((DATA / 'recovered_tables/index.json').read_text())
        for name, relative, cells, expected in CONFIG:
            table, checks = regenerate(name, relative, cells, expected)
            tables[name] = table
            processed.extend(checks)
            rel = 'working/figures/' + relative
            entries[rel] = {'file': f'data/genomics/generated_tables/{name}.json', 'complete': True,
                            'rows': len(table), 'cells': cells, 'notebook_sha256': sha(ROOT / rel)}
            if rel in recovered:
                reference = pd.read_json(ROOT / recovered[rel]['file'], orient='table')
                keys = ['Pop', 'abs_position' if 'abs_position' in reference else 'position']
                count = compare_table(reference, table, keys, displayed=True)
                comparisons.append({'dataset': name, 'saved_rows_matched': count, 'regenerated_rows': len(table)})
            print(f'Completed {name}: {len(table)} rows', flush=True)
        for relative, name in ALIASES.items():
            source = next('working/figures/' + r for n, r, _, _ in CONFIG if n == name)
            # The duplicate notebooks must share every executed data-processing cell.
            alias = json.loads((ROOT / 'working/figures' / relative).read_text())
            original = json.loads((ROOT / source).read_text())
            for cell in entries[source]['cells']:
                if alias['cells'][cell]['source'] != original['cells'][cell]['source']:
                    raise ValueError(f'Duplicate notebook data code differs: {relative}, cell {cell}')
            entries['working/figures/' + relative] = {**entries[source], 'executed_notebook': source,
                'notebook_sha256': sha(ROOT / 'working/figures' / relative)}
        if len(processed) != 17:
            raise ValueError(f'Expected 17 processed-export comparisons, found {len(processed)}')
        verify_inputs()  # Supplied inputs must remain byte-identical after the run.
        DEST.mkdir(exist_ok=True)
        for name, table in tables.items():
            table.to_json(DEST / f'{name}.json', orient='table', index=False, indent=2)
            csv_table = table.copy()
            csv_table['freq'] = csv_table['freq'].map(lambda x: json.dumps(vector(x)))
            csv_table.to_csv(DEST / f'{name}.csv', index=False)
        (DEST / 'index.json').write_text(json.dumps(entries, indent=2) + '\n')
        report = {'status': 'passed', 'mutation_input_files': 149, 'processed_exports_matched': processed,
                  'archive_sha256': manifest['archive_sha256'], 'saved_display_comparisons': comparisons,
                  'row_counts': {name: len(table) for name, table in tables.items()},
                  'comparison_rules': 'Mutation keys exact; numeric tolerance 1e-12; saved display prefixes and empty-field formatting normalized. Original notebook frequency rounding retained.',
                  'inputs_unchanged': True}
        (DEST / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
        status.write_text(json.dumps(report, indent=2) + '\n')
    except BaseException as error:
        status.write_text(json.dumps({'status': 'failed', 'error': str(error)}, indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
