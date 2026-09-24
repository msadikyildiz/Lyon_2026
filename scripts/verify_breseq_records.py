"""Import and verify archived references and records for every study breseq run."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import shlex
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/genomics'
RECORDS = BASE / 'run_records'
PREFIX = 'Breseq_reference_and_run_records/'
SOURCE = PREFIX + 'project/greencenter/Toprak_lab/shared/adam/adam_sequencing/'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sequences(data, genbank=False):
    text = data.decode()
    if genbank:
        seqs = [re.sub('[^a-zA-Z]', '', block.split('ORIGIN', 1)[1]).upper()
                for block in text.split('//') if 'ORIGIN' in block]
    else:
        seqs = [''.join(block.splitlines()[1:]).replace(' ', '').upper()
                for block in text.split('>')[1:]]
    if not seqs or any(not s or re.search('[^ACGTRYSWKMBDHVN]', s) for s in seqs):
        raise ValueError('Missing or invalid reference sequence')
    return Counter((len(s), digest(s.encode())) for s in seqs)


def study_runs():
    with (BASE / 'reference/metadata_complete.csv').open(encoding='utf-8-sig') as stream:
        return {f"{r['FolderDate']}/{r['source_file']}" for r in csv.DictReader(stream)}


def import_archive(archive):
    staged = {}
    files = []
    with ZipFile(archive) as z:
        names = [n for n in z.namelist() if not n.endswith('/')]
        observed = {n[len(SOURCE + 'out/'):].rsplit('/data/reference.fasta', 1)[0]
                    for n in names if n.startswith(SOURCE + 'out/') and n.endswith('/data/reference.fasta')}
        assert observed == study_runs(), 'Archive run identities differ from study metadata'
        assert len(names) == 600, 'Unexpected archive inventory'
        for name in sorted(names):
            data = z.read(name)
            if name.startswith(SOURCE + 'out/'):
                relative = name[len(SOURCE):]
                suffix = relative.rsplit('/', 2)[-2:]
                assert '/'.join(suffix) in {'data/reference.fasta', 'data/summary.json', 'output/summary.json', 'output/log.txt'}
                if name.endswith('reference.fasta'):
                    target = 'reference/fasta/' + digest(data) + '.fasta'
                else:
                    target = relative
            elif name.endswith(('.gb', '.gbk')):
                target = 'reference/genbank/' + Path(name).name
            elif name == PREFIX + 'README.txt':
                target = 'run_records/supplied_README.txt'
            else:
                raise ValueError('Unexpected archive member: ' + name)
            if target in staged:
                assert staged[target] == data, target
            staged[target] = data
            files.append({'archive_path': name, 'path': target, 'bytes': len(data), 'sha256': digest(data)})
    # Refuse to overwrite any differing existing input.
    for target, data in staged.items():
        p = BASE / target
        if p.exists() and p.read_bytes() != data:
            raise ValueError('Existing input differs: ' + str(p))
    for target, data in staged.items():
        p = BASE / target
        p.parent.mkdir(parents=True, exist_ok=True)
        if not p.exists():
            p.write_bytes(data)
    manifest = {'archive_name': archive.name, 'archive_sha256': digest(archive.read_bytes()),
                'source': 'Breseq_reference_and_run_records.zip',
                'storage': 'Original bytes retained; identical per-run FASTA files share one physical copy.',
                'files': files}
    RECORDS.mkdir(exist_ok=True)
    (RECORDS / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


def verify():
    manifest = json.loads((RECORDS / 'manifest.json').read_text())
    checked = {}
    fasta_by_run = {}
    for item in manifest['files']:
        p = BASE / item['path']
        if item['path'] not in checked:
            data = p.read_bytes()
            checked[item['path']] = (len(data), digest(data))
        assert checked[item['path']] == (item['bytes'], item['sha256']), str(p)
        if item['archive_path'].endswith('/data/reference.fasta'):
            run = item['archive_path'][len(SOURCE + 'out/'):].rsplit('/data/', 1)[0]
            fasta_by_run[run] = item['path']
    assert set(fasta_by_run) == study_runs()
    reference_cache = {}
    rows = []
    for run in sorted(study_runs()):
        output = BASE / 'out' / run / 'output'
        data = BASE / 'out' / run / 'data'
        assert (output / 'summary.json').read_bytes() == (data / 'summary.json').read_bytes(), run
        summary = json.loads((output / 'summary.json').read_text())
        log = (output / 'log.txt').read_text()
        commands = [s.strip() for s in log.splitlines() if s.strip().startswith('breseq ')]
        assert len(commands) == 1, run
        command = shlex.split(commands[0])
        assert '--polymorphism-prediction' in command
        assert command[command.index('--num-processors') + 1] == '16'
        original_ref = command[command.index('-r') + 1]
        reference = 'reference/genbank/' + Path(original_ref).name
        fasta = fasta_by_run[run]
        for path, gb in [(reference, True), (fasta, False)]:
            if path not in reference_cache:
                reference_cache[path] = sequences((BASE / path).read_bytes(), gb)
        assert reference_cache[reference] == reference_cache[fasta], run
        expected_length = sum(length * count for (length, _), count in reference_cache[fasta].items())
        assert summary['references']['total_length'] == expected_length, run
        options = summary['options']
        assert options['mutation_identification']['polymorphism_prediction'] is True
        assert options['workflow']['num_processors'] == 16
        assert options['mutation_identification']['polymorphism_frequency_cutoff'] == 0.05
        read_files = list(summary['reads']['read_file'])
        assert len(read_files) == 1 and read_files[0].endswith('_R1_001'), run
        rows.append({'run': run, 'recorded_date': log.splitlines()[0],
                     'recorded_reference_path': original_ref, 'archived_genbank': reference,
                     'reference_fasta': fasta, 'reference_length': expected_length,
                     'summary_json': (output / 'summary.json').relative_to(BASE).as_posix(),
                     'log': (output / 'log.txt').relative_to(BASE).as_posix(),
                     'recorded_command': commands[0]})
    with (RECORDS / 'runs.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {'study_runs': len(rows), 'archive_files_verified': len(manifest['files']),
              'unique_files_verified': len(checked), 'unique_run_fastas': len(set(fasta_by_run.values())),
              'reference_runs': dict(Counter(Path(r['archived_genbank']).name for r in rows)),
              'all_reference_sequences_match': True, 'both_summaries_identical_for_every_run': True,
              'runs_with_one_R1_labelled_input': len(rows),
              'polymorphism_frequency_cutoff': 0.05,
              'breseq_executable_version': None,
              'reference_provenance': 'Archived GenBank copies match per-run FASTA sequence content. Original /work paths are unavailable; original annotation-file byte identity is not established.'}
    (RECORDS / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--import-archive', type=Path)
    args = parser.parse_args()
    if args.import_archive:
        import_archive(args.import_archive)
    verify()
