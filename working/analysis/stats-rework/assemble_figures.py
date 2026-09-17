"""Compose the ten figures from checked sources and explicit panel geometry.

assembly_layout.json records the panel layout. assembly_assets.json records the
unchanged panel sources, conservative white-margin crops, hashes and panel identities.
Only composition changes here; numerical results and individual panels are inputs.
"""
from collections import Counter
from pathlib import Path
import hashlib
import io
import json
import math

import fitz
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'working/figures-assembled'


def source_bytes(asset):
    data = (ROOT / asset['source']).read_bytes()
    assert hashlib.sha256(data).hexdigest() == asset['sha256'], asset['source']
    return data


def ink_box(position, text, size, font='hebo'):
    return fitz.Rect(position[0], position[1] - 0.72 * size,
                     position[0] + fitz.get_text_length(text, fontname=font, fontsize=size),
                     position[1] + 0.02 * size)


def validate(layout, inventory):
    """Reject altered coverage, distortion and overlapping or out-of-page content."""
    expected = inventory['figures']
    assets = inventory['assets']
    assert Counter(f['name'] for f in layout['figures']) == Counter(expected.keys())
    for figure in layout['figures']:
        width, height = figure['width'], figure['height']
        assert width > 0 and height > 0 and height / width <= 7.5 / 6.45
        page = fitz.Rect(0, 0, width, height)
        identity = lambda p: (p['asset'], p['panel'], p['display_label'])
        assert Counter(map(identity, figure['panels'])) == Counter(map(identity, expected[figure['name']]))
        boxes = [('title', ink_box(layout['title_position'], figure['title'], layout['title_fontsize']))]
        for panel in figure['panels']:
            rect = fitz.Rect(panel['rect'])
            assert all(math.isfinite(v) for v in rect) and not rect.is_empty
            assert math.isclose(rect.width / rect.height, assets[panel['asset']]['content_aspect'], rel_tol=1e-4)
            boxes.append((panel['asset'], rect))
            if panel['display_label']:
                assert panel['label_position'] is not None
                boxes.append(('label ' + panel['display_label'], ink_box(
                    panel['label_position'], panel['display_label'], layout['panel_label_fontsize'])))
            else:
                assert panel['label_position'] is None
        for note in figure.get('notes', []):
            boxes.append(('note', ink_box(note['position'], note['text'], note['fontsize'], 'helv')))
        for i, (name, box) in enumerate(boxes):
            assert page.contains(box), (figure['name'], name, 'outside page')
            for other_name, other in boxes[:i]:
                overlap = box & other
                assert overlap.is_empty or overlap.width < 0.01 or overlap.height < 0.01, (
                    figure['name'], name, other_name, 'overlap')


def use_print_panels(layout, inventory):
    """Use size-specific vector panels while retaining the original source inventory."""
    manifest = json.loads((HERE / 'out/print-panels/manifest.json').read_text())
    assert manifest['layout_sha256'] == hashlib.sha256((HERE / 'assembly_layout.json').read_bytes()).hexdigest()
    for source, digest in manifest['inputs'].items():
        assert hashlib.sha256((ROOT / source).read_bytes()).hexdigest() == digest, source
    for figure in layout['figures']:
        for panel in figure['panels']:
            original = panel['asset']
            key = figure['name'] + '::' + original
            if key not in manifest['panels']:
                continue
            replacement = figure['name'] + '__' + original
            inventory['assets'][replacement] = dict(manifest['panels'][key],
                original_source=inventory['assets'][original]['source'])
            panel['asset'] = replacement
            for expected in inventory['figures'][figure['name']]:
                if expected['asset'] == original:
                    expected['asset'] = replacement


def main():
    layout = json.loads((HERE / 'assembly_layout.json').read_text())
    inventory = json.loads((HERE / 'assembly_assets.json').read_text())
    use_print_panels(layout, inventory)
    validate(layout, inventory)
    assets = inventory['assets']
    used = {p['asset'] for f in layout['figures'] for p in f['panels']}
    data = {key: source_bytes(assets[key]) for key in used}
    OUT.mkdir(exist_ok=True)
    manifest, checks = [], []
    for figure in layout['figures']:
        doc = fitz.open()
        page = doc.new_page(width=figure['width'], height=figure['height'])
        page.insert_text(layout['title_position'], figure['title'], fontsize=layout['title_fontsize'], fontname='hebo')
        expected_text = Counter(figure['title'].split())
        expected_vectors = expected_images = 0
        for panel in figure['panels']:
            asset = assets[panel['asset']]
            rect, crop = fitz.Rect(panel['rect']), fitz.Rect(asset['content_box'])
            if asset['source'].endswith('.pdf'):
                with fitz.open(stream=data[panel['asset']], filetype='pdf') as src:
                    assert len(src) == 1 and src[0].rect.contains(crop)
                    page.show_pdf_page(rect, src, 0, clip=crop, keep_proportion=True)
                    expected_text.update(w[4] for w in src[0].get_text('words'))
                    expected_vectors += len(src[0].get_drawings())
            else:
                with Image.open(io.BytesIO(data[panel['asset']])) as image:
                    assert image.size == (asset['native_width'], asset['native_height'])
                    assert fitz.Rect(0, 0, *image.size).contains(crop)
                    cropped = image.crop(tuple(int(round(v)) for v in crop))
                    stream = io.BytesIO()
                    cropped.save(stream, format='PNG')
                    page.insert_image(rect, stream=stream.getvalue(), keep_proportion=True)
                    expected_images += 1
            if panel['display_label']:
                page.insert_text(panel['label_position'], panel['display_label'],
                                 fontsize=layout['panel_label_fontsize'], fontname='hebo')
                expected_text.update([panel['display_label']])
            manifest.append(dict(figure=figure['title'], panel=panel['panel'],
                                 display_label=panel['display_label'], source=asset['source'],
                                 source_sha256=asset['sha256'], rectangle=panel['rect'],
                                 original_source=asset.get('original_source', asset['source']),
                                 source_content_box=asset['content_box'],
                                 label_position=panel['label_position'],
                                 page_size=[figure['width'], figure['height']]))
        for note in figure.get('notes', []):
            page.insert_text(note['position'], note['text'], fontsize=note['fontsize'], fontname='helv')
            expected_text.update(note['text'].split())
        actual_text = Counter(w[4] for w in page.get_text('words'))
        assert actual_text == expected_text, (figure['name'], expected_text - actual_text, actual_text - expected_text)
        assert len(page.get_drawings()) == expected_vectors
        assert len(page.get_images()) == expected_images
        doc.set_metadata({'title': figure['title'],
                          'subject': 'Lyon et al. (2026). Source data and analysis code: https://github.com/msadikyildiz/Lyon_2026'})
        destination = OUT / (figure['name'] + '.pdf')
        doc.save(destination, garbage=4, deflate=True)
        # Render the saved PDF so media-box rounding agrees with downstream viewers.
        with fitz.open(destination) as saved:
            saved[0].get_pixmap(dpi=300).save(OUT / (figure['name'] + '.png'))
            saved[0].get_pixmap(dpi=100).save(OUT / (figure['name'] + '_preview.png'))
        checks.append(dict(figure=figure['name'], panels=len(figure['panels']),
                           vector_paths=expected_vectors, raster_panels=expected_images,
                           source_words=sum(expected_text.values()), text_preserved=True))
        doc.close()
    (OUT / 'assembly_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    with fitz.open() as combined:
        ordered = sorted(layout['figures'], key=lambda f: (
            f['name'].startswith('Supplementary'), int(f['name'].split('_')[-1])))
        for figure in ordered:
            with fitz.open(OUT / (figure['name'] + '.pdf')) as source:
                combined.insert_pdf(source)
        combined.set_toc([[1, f['title'], i + 1] for i, f in enumerate(ordered)])
        combined.set_metadata({'title': 'Lyon et al. (2026): figures'})
        combined.save(OUT / 'All_figures.pdf', garbage=4, deflate=True)
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
