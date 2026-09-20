"""Compare regenerated tables, workbook cells and PNGs with a separate reference checkout."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from openpyxl import load_workbook
RTOL, ATOL = 1e-8, 1e-10
OUT = Path('working/analysis/stats-rework/out')

def compare(reference, rebuilt):
    if reference.resolve() == rebuilt.resolve():
        raise ValueError('Reference and rebuilt directories must differ')
    names = {p.name for p in (reference/OUT).glob('*.csv')}
    if not names or names != {p.name for p in (rebuilt/OUT).glob('*.csv')}:
        raise ValueError('CSV output file set differs or is empty')
    for name in sorted(names):
        a, b = pd.read_csv(reference/OUT/name), pd.read_csv(rebuilt/OUT/name)
        if name == 'fig5a_descriptive.csv':
            for x, y in zip(a['values'], b['values']):
                np.testing.assert_allclose(list(map(float,x.split(';'))), list(map(float,y.split(';'))), rtol=RTOL, atol=ATOL)
            a, b = a.drop(columns='values'), b.drop(columns='values')
        pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=RTOL,atol=ATOL,obj=name)
    pngs = {}
    for relative in [OUT/'figures-final', OUT/'print-panels', OUT/'trajectories', Path('working/figures-assembled')]:
        files = {p.name for p in (reference/relative).glob('*.png')}
        if not files or files != {p.name for p in (rebuilt/relative).glob('*.png')}:
            raise ValueError(f'PNG file set differs or is empty: {relative}')
        for name in sorted(files):
            a, b = reference/relative/name, rebuilt/relative/name
            if hashlib.sha256(a.read_bytes()).digest() != hashlib.sha256(b.read_bytes()).digest():
                raise ValueError(f'PNG differs: {relative/name}')
        pngs[str(relative)] = len(files)
    genomic = Path('data/genomics/generated_tables')
    genomic_files = {p.name for p in (reference / genomic).glob('*') if p.is_file()}
    if not genomic_files or genomic_files != {p.name for p in (rebuilt / genomic).glob('*') if p.is_file()}:
        raise ValueError('Regenerated genomic file set differs or is empty')
    for name in genomic_files:
        if (reference / genomic / name).read_bytes() != (rebuilt / genomic / name).read_bytes():
            raise ValueError(f'Regenerated genomic output differs: {name}')
    a=load_workbook(reference/'working/source-data/Source Data.xlsx',read_only=True,data_only=False)
    b=load_workbook(rebuilt/'working/source-data/Source Data.xlsx',read_only=True,data_only=False)
    if a.sheetnames != b.sheetnames:
        raise ValueError('Workbook sheets differ')
    numeric = cells = 0
    for name in a.sheetnames:
        if (a[name].max_row,a[name].max_column)!=(b[name].max_row,b[name].max_column):
            raise ValueError(f'Sheet dimensions differ: {name}')
        for ar,br in zip(a[name],b[name]):
            for x,y in zip(ar,br):
                if x.data_type != y.data_type:
                    raise ValueError(f'Cell type differs: {name}!{x.coordinate}')
                if isinstance(x.value,(int,float)):
                    if not isinstance(y.value,(int,float)) or not np.isclose(x.value,y.value,rtol=RTOL,atol=ATOL,equal_nan=True):
                        raise ValueError(f'Numeric cell differs: {name}!{x.coordinate}')
                    numeric += 1
                elif x.value != y.value:
                    # Doubling-time values are also stored as semicolon-separated numeric text.
                    try:
                        xv=list(map(float,x.value.split(';')));yv=list(map(float,y.value.split(';')))
                    except (AttributeError,ValueError):
                        raise ValueError(f'Text cell differs: {name}!{x.coordinate}') from None
                    np.testing.assert_allclose(xv,yv,rtol=RTOL,atol=ATOL)
                cells += 1
    a.close();b.close()
    return {'status':'passed','csv_tables':len(names),'pngs_identical':pngs,'workbook_cells_checked':cells,'genomic_files_identical':len(genomic_files),
            'workbook_numeric_cells':numeric,'rtol':RTOL,'atol':ATOL,
            'comparison':'CSV values, PNG bytes, workbook values and cell types; PDF/SVG/XLSX container timestamps ignored'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference',type=Path);parser.add_argument('rebuilt',type=Path)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    # Clear any previous success report before a comparison that may fail.
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps({'status':'running'})+'\n')
    try:
        report=compare(args.reference,args.rebuilt)
    except Exception as error:
        if args.report:args.report.write_text(json.dumps({'status':'failed','error':str(error)},indent=2)+'\n')
        raise
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
