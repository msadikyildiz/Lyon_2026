"""Rebuild the available analysis and figures; optionally refit in an isolated copy."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
SCRIPTS = Path('working/analysis/stats-rework')
DATASETS = ['fig2_paplpc', 'fig3_plac', 'supp3_pcr', 'supp4_atec', 'fig5_mutants', 'supp6_unt']
STEPS = ['validate.py', 'final_stats.py', 'mdk_source.py', 'supp11_mdk.py',
         'manifest.py', 'report_results.py', 'survival.py', 'plot_final.py', 'source_data.py',
         'print_panels.py', 'assemble_figures.py']


def run(root, args, log):
    env = os.environ.copy()
    env.setdefault('MPLCONFIGDIR', str(Path.home() / '.cache' / 'matplotlib'))
    env['MPLBACKEND'] = 'Agg'
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[name] = '1'
    print('Running ' + ' '.join(map(str, args)), flush=True)
    with log.open('w') as stream:
        result = subprocess.run([sys.executable, *map(str, args)], cwd=root,
                                env=env, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Exit {result.returncode}; see {log}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refit', action='store_true', help='Refit all six datasets in an isolated copy, then run the full pipeline there.')
    parser.add_argument('--jobs', type=int, default=2, help='Concurrent datasets during refitting (default: 2).')
    parser.add_argument('--output', type=Path, help='New directory for the isolated refit run; must not exist.')
    args = parser.parse_args()
    if args.jobs < 1 or (args.output and not args.refit):
        parser.error('--jobs must be positive; --output requires --refit')
    start = time.time()
    root = ROOT
    if args.refit:
        root = (args.output or ROOT / 'runs' / ('refit-' + time.strftime('%Y%m%d-%H%M%S'))).resolve()
        if root.exists():
            parser.error(f'Refit output already exists: {root}')
        if ROOT == root or ROOT.is_relative_to(root):
            parser.error('Refit output cannot contain the source checkout')
        if root.is_relative_to(ROOT) and not root.is_relative_to(ROOT / 'runs'):
            parser.error('Within this checkout, place refit outputs under runs/')
        shutil.copytree(ROOT, root, ignore=shutil.ignore_patterns('.git', 'runs', '__pycache__', '.DS_Store', '.venv'))
        # Remove only generated directories inside the newly created refit copy.
        for relative in [SCRIPTS / 'out', Path('working/source-data'), Path('working/figures-assembled'),
                         Path('data/genomics/generated_tables')]:
            if (root / relative).exists():
                shutil.rmtree(root / relative)
        cache = root / SCRIPTS / 'cache'
        cache.rename(root / 'reference-cache')
        cache.mkdir()
    logs = root / 'runs' / 'logs'
    logs.mkdir(parents=True, exist_ok=True)
    result_file = root / 'runs/reproduction.json'
    result_file.write_text(json.dumps({'status': 'running', 'started': time.strftime('%Y-%m-%dT%H:%M:%S%z')})+'\n')
    try:
        if args.refit:
            def refit(name):
                run(root, [SCRIPTS / 'build_all.py', name], logs / (name + '.log'))
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                list(pool.map(refit, DATASETS))
            run(root, ['scripts/compare_caches.py', str(root / 'reference-cache'), str(root / SCRIPTS / 'cache'), str(root / 'runs/refit-comparison.json')], logs / 'compare-refits.log')
        run(root, ['-m', 'unittest', 'discover', '-s', SCRIPTS / 'tests', '-v'], logs / 'tests.log')
        run(root, ['scripts/reproduce_genomics.py'], logs / 'genomics.log')
        for step in STEPS:
            run(root, [SCRIPTS / step], logs / (step + '.log'))
            if step == 'final_stats.py':
                run(root, ['scripts/check_primary.py'], logs / 'primary-regression.log')
        run(root, ['scripts/reproduce_trajectories.py'], logs / 'trajectories.log')
        if args.refit:
            run(root, ['scripts/compare_outputs.py', str(ROOT), str(root), '--report', str(root/'runs/output-comparison.json')], logs/'compare-outputs.log')
    except BaseException as error:
        result_file.write_text(json.dumps({'status': 'failed', 'error': str(error), 'seconds': round(time.time()-start, 2)}, indent=2)+'\n')
        raise
    report = {'mode': 'raw-refit' if args.refit else 'cached', 'python': platform.python_version(),
              'platform': platform.platform(), 'seconds': round(time.time()-start, 2),
              'steps': ['scripts/reproduce_genomics.py'] + STEPS + ['scripts/reproduce_trajectories.py'], 'status': 'passed',
              'coverage': 'docs/FIGURE_COVERAGE.md',
              'raw_inputs': {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted((root / SCRIPTS / 'biohpc-pull/data').glob('*.xlsx'))}}
    (root / 'runs/reproduction.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'Completed: {root / "working/figures-assembled/All_figures.pdf"}', flush=True)
    print(f'Validation: {root / "runs/reproduction.json"}', flush=True)

if __name__ == '__main__':
    main()
