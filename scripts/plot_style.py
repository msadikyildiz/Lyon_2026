"""Use the reference typeface when installed, otherwise Matplotlib's bundled serif."""
import os
from matplotlib import font_manager


def figure_font():
    requested = os.environ.get('LYON_PLOT_FONT', 'Times New Roman')
    try:
        font_manager.findfont(requested, fallback_to_default=False)
    except ValueError:
        return 'STIXGeneral'
    return requested
