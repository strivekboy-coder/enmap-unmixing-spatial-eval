"""Explicit raster alignment, reference rasterization and spectral response integration."""
from __future__ import annotations

import hashlib
import json
import logging
import re
import sqlite3
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from rasterio.features import rasterize
from shapely import make_valid
from shapely.validation import explain_validity

LOG = logging.getLogger(__name__)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def extract_data(cfg, root):
    directory = root / cfg['data']['directory']
    directory.mkdir(parents=True, exist_ok=True)
    required = ['enmap_berlin.tif', 'library_berlin.gpkg', 'landcover_berlin.gpkg']
    if not all((directory / name).exists() for name in required):
        with ZipFile(cfg['data']['archive']) as archive:
            for info in archive.infolist():
                if not info.is_dir():
                    # Flatten known tutorial inputs; never trust archive paths as destinations.
                    filename = Path(info.filename).name
                    if filename in required + ['hymap_berlin.tif']:
                        (directory / filename).write_bytes(archive.read(info))
    return directory


def unit_factor(unit):
    if str(unit).lower() in ('micrometers', 'micrometer', 'um', 'µm'):
        return 1000.0
    if str(unit).lower() in ('nanometers', 'nanometer', 'nm'):
        return 1.0
    raise ValueError(f'Unknown wavelength unit: {unit!r}')


def load_inputs(cfg, root, out):
    directory = extract_data(cfg, root)
    with rasterio.open(directory / 'enmap_berlin.tif') as src:
        image = src.read(masked=True).astype('float32').filled(np.nan)
        profile = src.profile.copy()
        wavelengths = np.array([
            float(src.tags(i)['wavelength']) * unit_factor(src.tags(i)['wavelength_units'])
            for i in range(1, src.count + 1)
        ])
        descriptions = list(src.descriptions)
    if not np.all(np.diff(wavelengths) > 0):
        raise ValueError('Image wavelengths must be strictly increasing.')
    classes = cfg['classes']
    spectra, labels, names = [], [], []
    matched = 0
    with sqlite3.connect(directory / 'library_berlin.gpkg') as con:
        for name, encoded, label in con.execute('SELECT name,profiles,level_1 FROM library_berlin ORDER BY fid'):
            obj = json.loads(encoded)
            axis = np.asarray(obj['x']) * unit_factor(obj['xUnit'])
            y = np.asarray(obj['y'], dtype=float)
            if len(axis) == len(wavelengths) and np.allclose(axis, wavelengths, atol=1e-4):
                matched += 1
            else:
                if wavelengths.min() < axis.min() or wavelengths.max() > axis.max():
                    raise ValueError(f'Library spectrum {name} does not cover image wavelengths.')
                y = np.interp(wavelengths, axis, y)
            spectra.append(y)
            labels.append(classes.index(label))
            names.append(name)
    library = np.asarray(spectra, dtype='float32')
    reasons = [[] for _ in wavelengths]
    spec = cfg['spectral']
    for i, wavelength in enumerate(wavelengths):
        if not spec['minimum_nm'] <= wavelength <= spec['maximum_nm']:
            reasons[i].append('configured spectral edge exclusion')
        for lo, hi in spec['exclude_intervals_nm']:
            if lo <= wavelength <= hi:
                reasons[i].append(f'water-vapour interval {lo}-{hi} nm')
        if not np.isfinite(library[:, i]).all():
            reasons[i].append('nonfinite library value')
        if not np.isfinite(image[i]).any():
            reasons[i].append('image band has no finite data')
    keep = np.array([not reason for reason in reasons])
    if keep.sum() < 10:
        raise ValueError('Too few usable spectral bands.')
    audit = pd.DataFrame({
        'band_1_based': np.arange(len(keep)) + 1, 'wavelength_nm': wavelengths,
        'kept': keep, 'reason': ['; '.join(reason) for reason in reasons],
        'source_description': descriptions,
    })
    audit.to_csv(out / 'band_audit.csv', index=False)
    # A nodata pixel is excluded from *both* sensor models, rather than dropping
    # an entire useful band because pixels outside the flight strip are missing.
    valid = np.isfinite(image[keep]).all(axis=0)
    raw_manifest = {
        'files': {name: {'sha256': sha256(directory / name), 'bytes': (directory / name).stat().st_size}
                  for name in ['enmap_berlin.tif', 'library_berlin.gpkg', 'landcover_berlin.gpkg']},
        'image_shape_bands_rows_cols': list(image.shape), 'crs': str(profile['crs']),
        'transform': list(profile['transform'])[:6], 'valid_image_pixels': int(valid.sum()),
        'library_size': len(library), 'class_counts': {c: int(np.sum(np.asarray(labels) == i)) for i, c in enumerate(classes)},
        'exact_wavelength_matches': matched, 'bands_retained': int(keep.sum()),
        'additional_bands_removed': int((~keep).sum()),
        'input_spectral_gaps_nm': [[float(a), float(b)] for a, b in zip(wavelengths[:-1], wavelengths[1:]) if b-a > 35],
    }
    (out / 'input_manifest.json').write_text(json.dumps(raw_manifest, indent=2), encoding='utf-8')
    LOG.info('Loaded %s, %d valid pixels, %d/%d bands, library=%s', image.shape, valid.sum(), keep.sum(), len(keep), raw_manifest['class_counts'])
    return image[keep], library[:, keep], np.array(labels), wavelengths[keep], valid, profile, raw_manifest


def reference_fractions(cfg, root, profile, valid_image, out):
    reference = gpd.read_file(root / cfg['data']['directory'] / 'landcover_berlin.gpkg').to_crs(profile['crs'])
    original = reference.copy()
    geometry_audit = []
    for index, row in reference.iterrows():
        before = row.geometry
        after = before
        if not before.is_valid:
            if not cfg['reference']['repair_invalid_geometry']:
                raise ValueError('Invalid reference geometry; repair is disabled.')
            after = make_valid(before)
            change = abs(after.area-before.area)/before.area
            if change > cfg['reference']['maximum_relative_geometry_area_change']:
                raise ValueError(f'Geometry repair changes reference area by {change:.4%}.')
            if after.geom_type not in ['Polygon', 'MultiPolygon']:
                raise ValueError(f'Unexpected repaired geometry type: {after.geom_type}')
            reference.loc[index, 'geometry'] = after
        geometry_audit.append({'feature': int(index), 'class': row.level_1,
                               'original_validity': explain_validity(before),
                               'area_before_m2': before.area, 'area_after_m2': after.area,
                               'relative_area_change': (after.area-before.area)/before.area})
    pd.DataFrame(geometry_audit).to_csv(out/'reference_geometry_audit.csv', index=False)
    ratio = abs(profile['transform'].a) / cfg['reference']['fine_resolution_m']
    if not np.isclose(ratio, round(ratio)):
        raise ValueError('Reference fine grid must divide the EnMAP grid exactly.')
    ratio = int(round(ratio))
    h, w = profile['height'], profile['width']
    classes = cfg['classes']
    shapes = [(row.geometry, classes.index(row.level_1) + 1) for row in reference.itertuples()]
    fine = rasterize(shapes, out_shape=(h*ratio, w*ratio),
                     transform=profile['transform'] * Affine.scale(1/ratio, 1/ratio),
                     fill=0, dtype='uint8', all_touched=cfg['reference']['all_touched'])
    original_fine = rasterize([(row.geometry, classes.index(row.level_1)+1) for row in original.itertuples()],
                              out_shape=fine.shape, transform=profile['transform']*Affine.scale(1/ratio, 1/ratio),
                              fill=0, dtype='uint8', all_touched=cfg['reference']['all_touched'])
    changed = int(np.sum(fine != original_fine))
    (out/'reference_geometry_raster_audit.json').write_text(json.dumps({'changed_3m_subpixels': changed,
                                                                     'total_subpixels': int(fine.size)}, indent=2), encoding='utf-8')
    LOG.info('Geometry repairs: %d invalid features; %d raster subpixels changed',
             int((~original.geometry.is_valid).sum()), changed)
    # Counts are divided by the complete 30 m pixel area, matching the official
    # algorithm; uncovered subpixels are not silently redistributed among classes.
    fractions = np.array([(fine == i+1).reshape(h, ratio, w, ratio).mean(axis=(1, 3))
                          for i in range(len(classes))], dtype='float32')
    coverage = fractions.sum(axis=0)
    valid = valid_image & (coverage >= cfg['reference']['minimum_polygon_coverage'] - 1e-7)
    if valid.sum() < 100:
        raise ValueError('Insufficient reference support.')
    write_raster(out / 'reference_fractions.tif', np.where(valid[None], fractions, np.nan), profile, classes)
    write_raster(out / 'polygon_coverage.tif', coverage[None], profile, ['polygon_coverage'])
    LOG.info('Reference: %d valid pixels at >= %.0f%% coverage', valid.sum(), 100*cfg['reference']['minimum_polygon_coverage'])
    return fractions, valid, reference


def read_srf_xlsx(path, platform, bands):
    """Read numeric SRF cells directly from the official OOXML workbook."""
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with ZipFile(path) as archive:
        strings = [''.join(node.itertext()) for node in ET.fromstring(archive.read('xl/sharedStrings.xml'))]
        workbook = ET.fromstring(archive.read('xl/workbook.xml'))
        sheets = workbook.find('s:sheets', ns)
        sheet = next(s for s in sheets if s.get('name') == f'Spectral Responses ({platform})')
        rid = sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
        relations = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        target = next(r.get('Target') for r in relations if r.get('Id') == rid)
        sheet_path = target.lstrip('/') if target.startswith('/') else 'xl/' + target
        rows = ET.fromstring(archive.read(sheet_path)).find('s:sheetData', ns)
        header, numeric = {}, []
        for row in rows:
            cells = {}
            for cell in row:
                node = cell.find('s:v', ns)
                if node is None:
                    continue
                col = re.match('[A-Z]+', cell.get('r')).group()
                cells[col] = strings[int(node.text)] if cell.get('t') == 's' else float(node.text)
            if not header:
                header = cells
                columns = [next(k for k, value in header.items() if value == f'{platform}_SR_AV_{band}') for band in bands]
            elif isinstance(cells.get('A'), (float, int)):
                values = [cells.get('A')] + [cells.get(col, np.nan) for col in columns]
                if np.isfinite(values).all():
                    numeric.append(values)
    values = np.asarray(numeric)
    if len(values) < 100:
        raise ValueError('Could not parse the official SRF workbook.')
    return values[:, 0], values[:, 1:]


def response_weights(wavelengths, response_axis, responses, max_gap):
    """Integrate a piecewise-linear spectrum, without bridging missing regions.

    Return normalized band weights and the proportion of full SRF area covered
    by the input spectrum. Flat spectra must remain flat after convolution.
    """
    n = len(wavelengths)
    right = np.searchsorted(wavelengths, response_axis, side='right')
    right = np.clip(right, 1, n-1)
    left = right - 1
    span = wavelengths[right] - wavelengths[left]
    supported = ((response_axis >= wavelengths[0]) & (response_axis <= wavelengths[-1]) & (span <= max_gap))
    t = (response_axis - wavelengths[left]) / span
    area = responses * np.gradient(response_axis)[:, None]
    weights = np.zeros((n, responses.shape[1]))
    np.add.at(weights, left[supported], area[supported] * (1-t[supported, None]))
    np.add.at(weights, right[supported], area[supported] * t[supported, None])
    coverage = weights.sum(axis=0) / area.sum(axis=0)
    weights /= weights.sum(axis=0)
    return weights.astype('float32'), coverage


def sensor_weights(cfg, root, wavelengths, out):
    spectral = cfg['spectral']
    path = root / spectral['srf_file']
    method = f"ESA SRF v4.0, {spectral['platform']}, linear interpolation and response integration"
    if not path.exists():
        from urllib.request import urlopen
        try:
            with urlopen(spectral['srf_url'], timeout=60) as stream:
                path.write_bytes(stream.read())
        except Exception as exc:
            if not spectral['allow_gaussian_fallback']:
                raise
            LOG.warning('Official SRF unavailable (%s); configured Gaussian fallback', exc)
    if path.exists():
        axis, responses = read_srf_xlsx(path, spectral['platform'], spectral['bands'])
        (out/'srf_provenance.json').write_text(json.dumps({'url': spectral['srf_url'], 'sha256': sha256(path),
                                                          'method': method}, indent=2), encoding='utf-8')
    else:
        # S2A approximate band centres/FWHM. This branch is explicitly labelled.
        pairs = [(492, 66), (560, 36), (665, 31), (704, 15), (740, 15), (783, 20), (833, 106), (865, 21), (1614, 91), (2202, 175)]
        axis = np.arange(300., 2601.)
        responses = np.array([np.exp(-.5*((axis-c)/(f/2.35482))**2) for c, f in pairs]).T
        method = 'Gaussian approximation using nominal S2A centre/FWHM'
    weights, coverage = response_weights(wavelengths, axis, responses, spectral['maximum_interpolation_gap_nm'])
    info = pd.DataFrame({'band': spectral['bands'], 'response_coverage': coverage,
                         'effective_wavelength_nm': wavelengths @ weights, 'method': method})
    info.to_csv(out / 's2_response_audit.csv', index=False)
    pd.DataFrame(weights, index=wavelengths, columns=spectral['bands']).to_csv(out / 's2_convolution_weights.csv', index_label='wavelength_nm')
    if np.min(coverage) < spectral['minimum_srf_coverage']:
        raise ValueError(f'Insufficient spectral coverage for S2: {info.to_dict("records")}')
    LOG.info('S2 convolution: %s; minimum response coverage %.3f', method, min(coverage))
    return weights, info


def write_raster(path, array, profile, names):
    metadata = profile.copy()
    metadata.update(driver='GTiff', count=array.shape[0], dtype='float32', nodata=-9999., compress='deflate')
    with rasterio.open(path, 'w', **metadata) as dst:
        dst.write(np.where(np.isfinite(array), array, -9999.).astype('float32'))
        dst.descriptions = tuple(names)
