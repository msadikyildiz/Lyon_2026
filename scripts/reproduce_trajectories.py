"""Execute the six original trajectory plotting notebooks with portable input/output paths."""
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', str(Path.home()/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ['working/figures/Figure 1/D-F/PAsurvival.ipynb', 'working/figures/Figure 1/D-F/PLsurvival.ipynb', 'working/figures/Figure 1/D-F/PCsurvival.ipynb', 'working/figures/Figure 3/A - survival/PLACsurvival.ipynb', 'working/figures/Supplemental Figure 4 - Pb/A - Survival/Pbsurvival.ipynb', 'working/figures/Supplemental Figure 3 - CefR/A - ResistanceEvo/pcr_evolution.ipynb']

def main():
    for rel in NOTEBOOKS:
        notebook = json.loads((ROOT/rel).read_text())
        namespace = {'__name__': '__main__'}
        for i, cell in enumerate(notebook['cells']):
            if cell['cell_type'] == 'code':
                exec(compile(''.join(cell['source']), f'{rel}:cell-{i}', 'exec'), namespace)
        plt.close('all')
        print('Completed ' + rel, flush=True)

if __name__ == '__main__':
    main()
