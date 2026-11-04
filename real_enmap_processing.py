"""Quality-controlled first look and frozen-library transfer to real EnMAP.

All numerical diagnostics, native subsets, models and GeoTIFFs remain in data/.
The only publishable artifacts are this code, two summary figures and one report.
"""
from __future__ import annotations

import json
import math
import sqlite3
from pathlib import Path

import geopandas as gpd
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
from rasterio.windows import Window, from_bounds

from data_io import write_raster
from real_enmap_io import StopWork, utc_now, write_json


def load_context(cfg, root):
    raw = root / cfg['raw_directory']
    meta = json.loads((raw / 'metadata_summary.json').read_text(encoding='utf-8'))
    with np.load(root / 'results/cache/prepared.npz') as source:
        prepared = {key: source[key] for key in source.files}
    with rasterio.open(root / 'data/enmap_berlin.tif') as src:
        profile = src.profile.copy()
        bounds = src.bounds
    return raw, meta, prepared, profile, bounds


def interpolation_matrix(source_axis, target_axis, maximum_gap, fwhm=None, boundary_max_nm=0):
    source_axis, target_axis = np.asarray(source_axis), np.asarray(target_axis)
    if not np.all(np.diff(source_axis) > 0):
        raise StopWork('Real sensor wavelengths are not strictly increasing; do not silently reorder detector bands.')
    weights = np.zeros((len(source_axis), len(target_axis)), dtype='float32')
    audit = []
    for column, wavelength in enumerate(target_axis):
        method = 'linear'
        right = int(np.searchsorted(source_axis, wavelength))
        if right < len(source_axis) and abs(source_axis[right] - wavelength) < 1e-6:
            left = right
        else:
            left = right - 1
            if left < 0 or right >= len(source_axis):
                raise StopWork('A Stage A wavelength is outside the real sensor support; no extrapolation is allowed.')
            if source_axis[right] - source_axis[left] > maximum_gap:
                nearest = min((left, right), key=lambda i: abs(source_axis[i] - wavelength))
                distance = abs(source_axis[nearest] - wavelength)
                if fwhm is None or distance > min(boundary_max_nm, fwhm[nearest] / 2):
                    raise StopWork('A Stage A wavelength crosses an unsupported real spectral gap.')
                left = right = nearest
                method = 'nearest_gap_boundary_within_half_FWHM_and_configured_limit'
        if left == right:
            weights[left, column] = 1
        else:
            fraction = (wavelength - source_axis[left]) / (source_axis[right] - source_axis[left])
            weights[left, column], weights[right, column] = 1 - fraction, fraction
        audit.append({'target_nm': float(wavelength), 'left_band_1_based': left + 1,
                      'right_band_1_based': right + 1, 'left_nm': float(source_axis[left]),
                      'right_nm': float(source_axis[right]), 'method': method})
    np.testing.assert_allclose(weights.sum(axis=0), 1, atol=1e-7)
    return weights, audit


def crop_window(src, bounds):
    floating = from_bounds(*bounds, transform=src.transform)
    left, top = max(0, math.floor(floating.col_off) - 2), max(0, math.floor(floating.row_off) - 2)
    right = min(src.width, math.ceil(floating.col_off + floating.width) + 2)
    bottom = min(src.height, math.ceil(floating.row_off + floating.height) + 2)
    return Window(left, top, right - left, bottom - top)


def warp(array, source_profile, target_profile, method=Resampling.nearest, fill=0):
    destination = np.full((target_profile['height'], target_profile['width']), fill, dtype='float32')
    reproject(array.astype('float32'), destination,
              src_transform=source_profile['transform'], src_crs=source_profile['crs'],
              dst_transform=target_profile['transform'], dst_crs=target_profile['crs'],
              resampling=method, src_nodata=None, dst_nodata=-9999 if np.isfinite(fill) else np.nan, init_dest_nodata=False,
              num_threads=2)
    return destination


def quality_stage(cfg, root, deadline):
    raw, meta, prepared, target, bounds = load_context(cfg, root)
    source_axis = np.array([band['wavelength_nm'] for band in meta['bands']])
    weights, audit = interpolation_matrix(source_axis, prepared['wavelengths'],
                                         cfg['spectral_alignment']['maximum_source_gap_nm'],
                                         np.array([b['fwhm_nm'] for b in meta['bands']]),
                                         cfg['spectral_alignment']['gap_boundary_nearest_max_nm'])
    needed = np.flatnonzero(weights.sum(axis=1) > 0)
    with rasterio.open(raw / 'image.tif') as src:
        if src.count != len(source_axis) or src.nodata != meta['background_raw_dn']:
            raise StopWork('Image band count/nodata disagrees with downloaded XML.')
        window = crop_window(src, bounds)
        native_profile = src.profile.copy()
        native_profile.update(width=int(window.width), height=int(window.height),
                              transform=src.window_transform(window))
        counts = src.read(window=window).astype('float32')
        raster_scales = list(src.scales)
    gains = np.array([band['gain'] for band in meta['bands']], dtype='float32')
    offsets = np.array([band['offset'] for band in meta['bands']], dtype='float32')
    finite_band = counts != meta['background_raw_dn']
    reflectance = counts * gains[:, None, None] + offsets[:, None, None]
    reflectance[~finite_band] = np.nan
    layers = {}
    for key in ['quality_classes', 'quality_cloud', 'quality_cloud_shadow', 'quality_haze',
                'quality_cirrus', 'quality_snow', 'quality_testflags']:
        deadline.check()
        with rasterio.open(raw / (key + '.tif')) as src:
            if src.transform != target_from_native(native_profile, window) or src.crs != native_profile['crs']:
                raise StopWork('Quality and spectral source grids differ; explicit alignment is required.')
            layers[key] = src.read(1, window=window)
    with rasterio.open(raw / 'defective_pixel_mask.tif') as src:
        if src.count != len(source_axis):
            raise StopWork('Defective-pixel mask band count differs from XML.')
        defect_cube = src.read(window=window) != 0
        defective = np.any(defect_cube[needed], axis=0)
    inside = np.isin(layers['quality_classes'], [1, 2])
    finite = np.all(finite_band[needed], axis=0)
    conditions = {key.removeprefix('quality_'): layers[key] != 0 for key in
                  ['quality_cloud', 'quality_cloud_shadow', 'quality_haze', 'quality_cirrus', 'quality_snow']}
    conditions['defective_in_model_support'] = defective
    conditions['any_quality_test_flag'] = layers['quality_testflags'] != 0
    atmospheric = inside.copy()
    for name in ('cloud', 'cloud_shadow', 'haze', 'cirrus', 'snow'):
        if cfg['quality']['reject_' + ('any_cirrus' if name == 'cirrus' else name)]:
            atmospheric &= ~conditions[name]
    if cfg['quality'].get('require_testflags_byte_zero', True):
        atmospheric &= ~conditions['any_quality_test_flag']
    elif cfg['quality']['require_nominal_testflags']:
        atmospheric &= (layers['quality_testflags'] & 3) == 0
    clear = atmospheric & finite
    if cfg['quality']['reject_defective_in_model_bands']:
        clear &= ~defective
    diagnostic_reflectance = reflectance.copy()
    diagnostic_reflectance[defect_cube | ~atmospheric[None]] = np.nan
    nearest_inside = warp(inside, native_profile, target) > .5
    native_clear_on_target = warp(clear, native_profile, target) > .5
    atmospheric_on_target = warp(atmospheric, native_profile, target) > .5
    # Accept a bilinear output only when every contributing source pixel passed QA.
    support_weight = warp(clear, native_profile, target, Resampling.bilinear)
    usable = support_weight >= 1 - 1e-6
    data_on_target = np.full((len(needed), target['height'], target['width']), np.nan, dtype='float32')
    for index, band in enumerate(needed):
        deadline.check()
        source = np.where(clear, reflectance[band], np.nan).astype('float32')
        reproject(source, data_on_target[index], src_transform=native_profile['transform'],
                  src_crs=native_profile['crs'], dst_transform=target['transform'], dst_crs=target['crs'],
                  src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.bilinear, num_threads=2)
    aligned = np.einsum('bhw,bt->thw', data_on_target, weights[needed], optimize=True)
    usable &= np.isfinite(aligned).all(axis=0)
    aligned[:, ~usable] = np.nan
    fractions = {}
    denominator = target['height'] * target['width']
    for name, condition in conditions.items():
        on_grid = warp(condition & inside, native_profile, target) > .5
        fractions[name] = {'pixels': int(on_grid.sum()), 'share_of_full_tutorial_grid': float(on_grid.mean())}
    inside_fraction = float(nearest_inside.mean())
    flags = layers['quality_testflags']
    flags_counts = {str(value): int(count) for value, count in zip(*np.unique(flags[inside], return_counts=True))}
    stats = {'time_utc': utc_now(), 'full_tutorial_grid_pixels': denominator,
             'scene_footprint_pixels_on_tutorial_grid': int(nearest_inside.sum()),
             'scene_footprint_share': inside_fraction,
             'quality_clear_nearest_pixels': int(native_clear_on_target.sum()),
             'atmospheric_and_testflag_clear_nearest_pixels': int(atmospheric_on_target.sum()),
             'atmospheric_and_testflag_clear_share': float(atmospheric_on_target.mean()),
             'conservative_usable_pixels': int(usable.sum()),
             'usable_share_of_full_tutorial_grid': float(usable.mean()),
             'tutorial_valid_image_pixels': int(prepared['valid_image'].sum()),
             'usable_pixels_on_tutorial_valid_image': int((usable & prepared['valid_image']).sum()),
             'usable_share_of_tutorial_valid_image': float((usable & prepared['valid_image']).sum() / prepared['valid_image'].sum()),
             'usable_reference_pixels': int((usable & prepared['valid']).sum()),
             'individual_flags_overlap_and_are_not_additive': fractions,
             'source_quality_testflag_values': flags_counts,
             'quality_rule': 'Class land/water; no cloud, shadow, haze, any cirrus, snow, defect in spectral support, or nonzero quality-test flag; all bilinear contributors valid.',
             'real_source_bands_needed_for_177_features': len(needed),
             'real_source_band_count': len(source_axis),
             'image_raster_scale_values': sorted(set(raster_scales)),
             'reflectance_scaling_applied': 'XML gain/offset exactly once to unscaled Rasterio DN',
             'spatial_resampling': 'bilinear to the exact tutorial grid; conservative QA support threshold 1-1e-6',
             'spectral_resampling': 'piecewise linear to frozen 177-band tutorial centres; no extrapolation or gap >35nm',
             'gap_boundary_nearest_mappings': [row for row in audit if row['method'] != 'linear'],
             'no_reference_tuning_or_sum_to_one': True}
    write_json(raw / 'quality_summary.json', stats)
    write_json(raw / 'spectral_alignment.json', audit)
    pd.DataFrame(audit).to_csv(raw / 'spectral_alignment.csv', index=False)
    np.savez_compressed(raw / 'quality_prepared.npz', aligned_reflectance=aligned, usable=usable,
                        source_reflectance=reflectance, source_clear=clear, source_needed=needed,
                        source_diagnostic_reflectance=diagnostic_reflectance, source_atmospheric_clear=atmospheric,
                        atmospheric_on_target=atmospheric_on_target,
                        source_axis_nm=source_axis, weights=weights, support_weight=support_weight,
                        native_transform=np.array(list(native_profile['transform'])[:6]))
    write_raster(raw / 'usable_mask.tif', usable[None].astype('float32'), target, ['conservative_usable'])
    print(json.dumps({k: v for k, v in stats.items() if k not in ['individual_flags_overlap_and_are_not_additive']}, indent=2), flush=True)
    refresh_report(cfg, root)
    return stats


def target_from_native(native_profile, window):
    return native_profile['transform'] * Affine.translation(-window.col_off, -window.row_off)


def historical_purity(root, profile, labels):
    polygons = gpd.read_file(root / 'data/landcover_berlin.gpkg').to_crs(profile['crs'])
    shapes = [(row.geometry, labels.index(row.level_3) + 1) for row in polygons.itertuples()
              if row.level_3 in labels]
    fine = rasterize(shapes, out_shape=(profile['height'] * 10, profile['width'] * 10),
                     transform=profile['transform'] * Affine.scale(.1, .1), fill=0, dtype='uint8')
    return np.array([(fine == index + 1).reshape(profile['height'], 10, profile['width'], 10).mean(axis=(1, 3))
                     for index in range(len(labels))])


def median_available(array):
    result = np.full(array.shape[0], np.nan)
    present = np.isfinite(array).any(axis=1)
    result[present] = np.nanmedian(array[present], axis=1)
    return result


def plot_supported(ax, wavelength, spectrum, **options):
    cuts = np.r_[0, np.flatnonzero(np.diff(wavelength) > 35) + 1, len(wavelength)]
    for start, stop in zip(cuts[:-1], cuts[1:]):
        ax.plot(wavelength[start:stop], spectrum[start:stop], **options)
        options.pop('label', None)


def spectra_stage(cfg, root, deadline):
    raw, meta, prepared, profile, bounds = load_context(cfg, root)
    with np.load(raw / 'quality_prepared.npz') as source:
        aligned, usable = source['aligned_reflectance'], source['usable']
        native, native_axis, native_clear = source['source_diagnostic_reflectance'], source['source_axis_nm'], source['source_atmospheric_clear']
        atmospheric_on_target = source['atmospheric_on_target']
        native_transform = Affine(*source['native_transform'])
    labels = cfg['spectra']['labels']
    purity = historical_purity(root, profile, labels)
    library = []
    with sqlite3.connect(root / 'data/library_berlin.gpkg') as connection:
        for name, encoded, level3 in connection.execute('SELECT name,profiles,level_3 FROM library_berlin ORDER BY fid'):
            obj = json.loads(encoded)
            library.append((name, level3, np.asarray(obj['y']) / 10000))
    wave = prepared['wavelengths']
    # A spectral roughness proxy is diagnostic, never used to tune or replace the RF.
    median = median_available(native[:, native_clear])
    rough = np.full(len(native_axis), np.nan)
    for band in range(1, len(native_axis) - 1):
        if max(native_axis[band] - native_axis[band-1], native_axis[band+1] - native_axis[band]) > 35:
            continue
        expected = np.interp(native_axis[band], native_axis[[band-1, band+1]], median[[band-1, band+1]])
        rough[band] = abs(median[band] - expected) / max(abs(median[band]), .02)
    finite_share = np.isfinite(native[:, native_clear]).mean(axis=1)
    noisy = (rough > cfg['spectra']['diagnostic_second_difference_relative_threshold']) | (finite_share < cfg['spectra']['diagnostic_minimum_band_finite_fraction'])
    pd.DataFrame({'wavelength_nm': native_axis, 'relative_second_difference_proxy': rough,
                  'finite_fraction_on_clear_pixels': finite_share, 'diagnostic_flag': noisy}).to_csv(raw / 'noise_diagnostics.csv', index=False)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    examples, counts, signature_diagnostics = [], {}, {}
    colors = ['#555555', '#2d8c4b', '#b78245', '#2279b5']
    red, green, blue = [int(np.argmin(abs(native_axis - value))) for value in (665, 560, 490)]
    native_profile = profile.copy()
    native_profile.update(transform=native_transform, height=native.shape[1], width=native.shape[2])
    rgb = np.stack([warp(native[b], native_profile, profile, fill=np.nan) for b in (red, green, blue)], axis=-1)
    rgb = np.clip(rgb / .25, 0, 1)
    rgb[~np.isfinite(rgb).all(axis=2)] = .93
    for index, label in enumerate(labels):
        deadline.check()
        candidates = np.argwhere((purity[index] >= cfg['spectra']['historical_purity']) & atmospheric_on_target & prepared['valid'])
        counts[label] = len(candidates)
        ax = axes.flat[index]
        members = [row[2] for row in library if row[1] == label]
        for spectrum in members:
            plot_supported(ax, wave, spectrum, color='#d0d0d0', lw=.6)
        if len(members):
            plot_supported(ax, wave, np.median(members, axis=0), color='#555555', ls=':', label='2009 library median')
        if len(candidates):
            # Fixed coordinate-stratified selection, without using prediction errors.
            order = np.linspace(0, len(candidates) - 1, min(cfg['spectra']['pixels_per_class'], len(candidates)), dtype=int)
            selected = candidates[order]
            real = []
            for row, column in selected:
                x, y = profile['transform'] * (int(column) + .5, int(row) + .5)
                column_native, row_native = ~native_transform * (x, y)
                real.append(native[:, int(math.floor(row_native)), int(math.floor(column_native))])
            real = np.asarray(real).T
            simulated = prepared['image'][:, selected[:, 0], selected[:, 1]] / 10000
            for row, column in selected:
                x, y = profile['transform'] * (int(column) + .5, int(row) + .5)
                examples.append({'historical_label_2009': label, 'row': int(row), 'column': int(column), 'x': x, 'y': y})
                axes[1, 1].plot(column, row, 'o', mfc='none', mec=colors[index], ms=6)
            actual_median = median_available(real)
            plot_supported(ax, native_axis, actual_median, color=colors[index], lw=1.5, label='2025 native bands; bad values masked')
            plot_supported(ax, wave, np.median(simulated, axis=1), color='#b44d4d', lw=1, ls='--', label='2009 simulated, same positions')
            r_index, n_index = [int(np.argmin(abs(native_axis-v))) for v in (665,850)]
            signature_diagnostics[label] = {'red_reflectance': float(actual_median[r_index]),
                                            'nir_reflectance': float(actual_median[n_index]),
                                            'ndvi_of_median_spectrum': float((actual_median[n_index]-actual_median[r_index]) /
                                                (actual_median[n_index]+actual_median[r_index]))}
        else:
            ax.text(.5, .55, 'No QA-clear pixels at >=95%\nhistorical purity; no spectrum invented',
                    ha='center', transform=ax.transAxes, fontsize=9)
        for lo, hi in cfg['spectra']['water_vapour_intervals_nm']:
            ax.axvspan(lo, hi, color='#94b5d9', alpha=.25)
        for wavelength in native_axis[noisy]:
            ax.axvline(wavelength, color='#d28534', alpha=.12, lw=.6)
        ax.set(title=f'{label}: historical pure support n={len(candidates)}', xlabel='Wavelength (nm)', ylabel='Reflectance')
        ax.set_xlim(420, 2450)
        ax.legend(fontsize=7, loc='upper right')
    axes[1, 1].imshow(rgb)
    axes[1, 1].set(title='2025 RGB, per-band QA; sampled historical labels', xlabel='Column', ylabel='Row')
    quality = json.loads((raw / 'quality_summary.json').read_text(encoding='utf-8'))
    entries = quality['individual_flags_overlap_and_are_not_additive']
    names = {'any_quality_test_flag': 'Any quality-test flag', 'defective_in_model_support': 'Defect in any model-support band',
             'cloud_shadow': 'Cloud shadow', 'cloud':'Cloud', 'cirrus':'Cirrus', 'haze':'Haze', 'snow':'Snow'}
    percentages = [entries[name]['share_of_full_tutorial_grid'] * 100 for name in entries]
    bars = axes[1, 2].barh([names[name] for name in entries], percentages, color='#667d94')
    axes[1, 2].bar_label(bars, labels=[f'{v:.2f}%' for v in percentages], padding=3, fontsize=8)
    axes[1, 2].set_xlim(0,100)
    axes[1, 2].set(xlabel='% of full tutorial grid (overlapping flags)', title=f'Complete-feature inference support: {usable.mean():.1%}')
    fig.suptitle('Real EnMAP first look: 2025 vs 2009 — 16-year gap; historical classes are not 2025 truth\n'
                 'Blue: water-vapour windows; orange: roughness/finite-support diagnostic, not measured sensor SNR', fontsize=12)
    figures = root / cfg['figure_directory']; figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / 'spectra_quality.png', dpi=170)
    plt.close(fig)
    pd.DataFrame(examples).to_csv(raw / 'spectrum_example_locations.csv', index=False)
    result = {'time_utc': utc_now(), 'historical_pure_clear_support': counts,
              'example_pixels': examples, 'diagnostic_flagged_source_bands': int(noisy.sum()),
              'median_spectrum_signature_diagnostics': signature_diagnostics,
              'noise_definition': 'Relative adjacent-band residual >0.1 or <95% finite support; diagnostic, not measured SNR and not an RF exclusion.',
              'class_interpretation': '2009 level_3 polygon labels, >=95% purity, coordinate-stratified samples on atmospheric/testflag-clear support; bad values masked per band. Not independently validated 2025 material labels.',
              'plot_interpretation': '2025 spectrum is measured reflectance at historical positions; no claim that old bare soil or roof remains unchanged.'}
    write_json(raw / 'spectra_summary.json', result)
    refresh_report(cfg, root)
    print(json.dumps(result, indent=2), flush=True)
    return result


def transfer_stage(cfg, root, deadline):
    raw, meta, prepared, profile, bounds = load_context(cfg, root)
    with np.load(raw / 'quality_prepared.npz') as source:
        reflectance, usable = source['aligned_reflectance'], source['usable']
    if not usable.any():
        raise StopWork('No conservatively usable real pixels; transfer is not meaningful.')
    bundle = joblib.load(raw / 'stage_a_models.joblib')
    replay = json.loads((raw / 'model_replay_audit.json').read_text(encoding='utf-8'))
    if replay['enmap_maximum_absolute_prediction_difference'] > 1e-7:
        raise StopWork('EnMAP estimator replay was not equivalent to Stage A.')
    features = np.ascontiguousarray(reflectance[:, usable].T * 10000, dtype='float32')
    spectral = {'EnMAP': features, 'S2-like': features @ bundle['s2_weights']}
    means, spreads, members = {}, {}, {}
    for sensor, inputs in spectral.items():
        values = np.empty((3, 4, len(inputs)), dtype='float32')
        for run in range(3):
            for target in range(4):
                deadline.check(5)
                values[run, target] = bundle['models'][(sensor, run, target)].predict(inputs)
        mean = np.full((4, *usable.shape), np.nan, dtype='float32')
        std = mean.copy()
        mean[:, usable], std[:, usable] = values.mean(axis=0), values.std(axis=0, ddof=0)
        means[sensor], spreads[sensor], members[sensor] = mean, std, values
        name = sensor.lower().replace('-', '_')
        write_raster(raw / f'{name}_real_fractions.tif', mean, profile, bundle['stage_a_config']['classes'])
        write_raster(raw / f'{name}_real_uncertainty.tif', std, profile, bundle['stage_a_config']['classes'])
    common = usable & prepared['valid']
    if not common.any():
        raise StopWork('No common QA-clear/reference pixels; maps exist but no historical error is claimed.')
    tables = []
    simulations = {}
    with np.load(root / 'results/cache/ensemble_predictions.npz') as saved:
        for sensor in spectral:
            grid = np.full((4, *common.shape), np.nan, dtype='float32')
            grid[:, prepared['valid_image']] = saved[sensor][:3].mean(axis=0)
            simulations[sensor] = grid
            for target, label in list(enumerate(bundle['stage_a_config']['classes'])) + [(None, 'macro')]:
                truth = prepared['reference'][:, common]
                real = means[sensor][:, common]
                simulated = grid[:, common]
                if target is not None:
                    truth, real, simulated = truth[target], real[target], simulated[target]
                tables.append({'sensor': sensor, 'class': label, 'n_common_reference_pixels': int(common.sum()),
                               'real_2025_vs_reference_2009_mae': float(np.abs(real - truth).mean()),
                               'simulated_2009_vs_reference_2009_mae_same_pixels': float(np.abs(simulated - truth).mean()),
                               'real_minus_simulated_mae': float(np.abs(real-truth).mean() - np.abs(simulated-truth).mean()),
                               'real_bias': float((real - truth).mean()),
                               'real_rmse': float(np.sqrt(((real-truth)**2).mean())),
                               'reference_positive_pixels': int((truth > 0).sum()) if target is not None else None,
                               'calendar_year_gap': 16})
    metrics = pd.DataFrame(tables)
    metrics.to_csv(raw / 'historical_comparison.csv', index=False)
    np.savez_compressed(raw / 'real_predictions.npz', usable=usable, common=common,
                        enmap_mean=means['EnMAP'], enmap_std=spreads['EnMAP'],
                        s2_like_mean=means['S2-like'], s2_like_std=spreads['S2-like'])
    fig, axes = plt.subplots(3, 4, figsize=(13, 10), constrained_layout=True)
    for target, label in enumerate(bundle['stage_a_config']['classes']):
        for row, image, title, limit in [(0, means['EnMAP'][target], 'EnMAP fractions', (0, 1)),
                                         (1, spreads['EnMAP'][target], 'Ensemble SD', (0, .15))]:
            ax = axes[row, target]
            draw = ax.imshow(image, vmin=limit[0], vmax=limit[1], cmap='viridis' if row != 2 else 'magma')
            ax.set_title(f'{label}: {title}', fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(draw, ax=ax, shrink=.7)
        measured = metrics[(metrics.sensor == 'EnMAP') & (metrics['class'] == label)].iloc[0]
        ax = axes[2, target]
        bars = ax.bar(['Real 2025', 'Simulated 2009'],
                      [measured.real_2025_vs_reference_2009_mae, measured.simulated_2009_vs_reference_2009_mae_same_pixels],
                      color=['#417b9a', '#b8b8b8'])
        ax.bar_label(bars, fmt='%.3f', fontsize=8, padding=3)
        ax.set_ylim(0, .4)
        ax.set_title(f'{label}: historical MAE on same 31 pixels', fontsize=9)
        ax.set_ylabel('MAE (fraction)')
        ax.set_xlabel(f'Reference-positive support: {int(measured.reference_positive_pixels)}')
    fig.suptitle(f'Frozen-library RF on real EnMAP 2025 — {common.sum()} common historical pixels; 16-year gap\n'
                 'Maps are local transfer diagnostics; historical differences include real land-cover change', fontsize=12)
    figures = root / cfg['figure_directory']; figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / 'transfer_summary.png', dpi=170)
    plt.close(fig)
    library = prepared['library']
    outside = (features < library.min(axis=0)) | (features > library.max(axis=0))
    result = {'time_utc': utc_now(), 'usable_inference_pixels': int(usable.sum()),
              'common_historical_reference_pixels': int(common.sum()), 'calendar_year_gap': 16,
              'historical_dominant_support': {name:int((prepared['reference'][:,common].argmax(axis=0)==i).sum())
                    for i,name in enumerate(bundle['stage_a_config']['classes'])},
              'feature_entries_outside_library_range_fraction': float(outside.mean()),
              'pixels_with_any_feature_outside_library_range_fraction': float(outside.any(axis=1).mean()),
              'real_reflectance_negative_entry_fraction': float((features < 0).mean()),
              'real_reflectance_above_one_entry_fraction': float((features > 10000).mean()),
              'fixed_model': 'Stage A library/synthetic RF; 3 members; same seeds/hyperparameters; no fitting on real/reference pixels',
              's2_like_branch': 'Same real EnMAP spectra convolved with Stage A weights, not an actual Sentinel-2 acquisition; sensitivity-only due to one nonidentical replay.',
              'geometry': 'Nominal georeferencing only; AROSICS has not been run.',
              'metrics': tables}
    write_json(raw / 'transfer_summary.json', result)
    refresh_report(cfg, root)
    print(metrics.to_string(index=False), flush=True)
    return result


def refresh_report(cfg, root, include_cross_mission=True):
    raw = root / cfg['raw_directory']
    def read(name):
        path = raw / name
        return json.loads(path.read_text(encoding='utf-8')) if path.exists() else None
    meta, quality, spectra, transfer = [read(name) for name in
        ('metadata_summary.json', 'quality_summary.json', 'spectra_summary.json', 'transfer_summary.json')]
    geometry = read('arosics_result.json')
    cross = read('cross_mission_summary.json') if include_cross_mission else None
    angle = meta['angles_deg_by_corner']
    quality_note = 'Quality calculation pending.'
    if quality:
        quality_note = (f"Atmospheric/testflag-clear nearest-grid share: {quality['atmospheric_and_testflag_clear_share']:.2%}. "
                        f"**Conservative complete-feature inference support: {quality['conservative_usable_pixels']:,}/{quality['full_tutorial_grid_pixels']:,} pixels "
                        f"({quality['usable_share_of_full_tutorial_grid']:.2%} of the full tutorial rectangle; "
                        f"{quality['usable_share_of_tutorial_valid_image']:.2%} of valid simulated-image support).** "
                        f"Historical reference overlap: {quality['usable_reference_pixels']} pixels. "
                        'Reject cloud, cloud shadow, haze, any cirrus, snow, defective source bands and any quality-test flag; '
                        'every bilinear contributor must pass. Individual flag shares overlap.')
    spectrum_note = 'Spectrum comparison pending.'
    if spectra:
        spectrum_note = ('Historical >=95%-pure atmospheric-clear support: ' + ', '.join(f'{k} {v}' for k,v in spectra['historical_pure_clear_support'].items())
                         + '. Samples use 2009 roof/low-vegetation/soil/water labels, not independently verified 2025 identities. '
                         'Water-vapour windows and a fixed spectral-roughness/finite-support proxy are marked; this proxy is not measured sensor SNR.')
        spectrum_note += (' The historical soil and water groups have 2025 median-spectrum NDVI '
                          f"{spectra['median_spectrum_signature_diagnostics']['soil']['ndvi_of_median_spectrum']:.2f} and "
                          f"{spectra['median_spectrum_signature_diagnostics']['water']['ndvi_of_median_spectrum']:.2f}; "
                          'their signatures warn against treating the old labels as current pure materials.')
    metric_note = 'Frozen RF inference and historical comparison pending.'
    if transfer:
        lines = ['| Fraction | Real 2025 vs 2009 MAE | Simulated 2009 MAE, same pixels |', '| --- | ---: | ---: |']
        for row in transfer['metrics']:
            if row['sensor'] == 'EnMAP':
                lines.append(f"| {row['class']} | {row['real_2025_vs_reference_2009_mae']:.4f} | {row['simulated_2009_vs_reference_2009_mae_same_pixels']:.4f} |")
        metric_note = (f"Frozen 3-member library RF applied to {transfer['usable_inference_pixels']:,} real QA-clear pixels. "
                       f"MAE below uses the identical {transfer['common_historical_reference_pixels']} reference pixels for both dates.\n\n"
                       + '\n'.join(lines))
        metric_note += ('\n\nDominant-reference support (impervious/vegetation/soil/water): '
                        + '/'.join(str(transfer['historical_dominant_support'][name]) for name in ['impervious','vegetation','soil','water'])
                        + '; these 31 spatially sparse pixels cannot establish generalization.')
    optional_note = 'Optional AROSICS/actual-S2 checks are not completed.'
    if geometry and geometry.get('input_mask_rule','').startswith('Independent'):
        reliable = [r for r in geometry['windows'] if r['accepted']]
        optional_note = (f"**Actual-S2 check:** Sentinel-2C 2025-05-12 (+1 day), scene cloud 8.00%. "
                         f"AROSICS {geometry['arosics_version']} accepted {len(reliable)} independent-mask 64x64 window "
                         f"(reliability {reliable[0]['reliability_percent']:.2f}%, an internal score). "
                         f"Correction to add to EnMAP origin: east {geometry['median_x_correction_m']:+.2f} m, "
                         f"north {geometry['median_y_correction_m']:+.2f} m. This is a local estimate, not scene-wide accuracy; RF results remain unshifted.")
        optional_note += (' Original ESA XML and red/NIR/SWIR JP2 comparisons (59,755 pixels/band) establish that '
                          'S2 COG DN already includes -1000 DN: effective reflectance offset is 0, avoiding double subtraction. No EnMAP-based radiometric fitting.')
    if cross:
        primary=[r for r in cross['metrics'] if r['geometry']=='arosics']
        optional_note += (f" Official S2C SRFs give 10 band comparisons, {min(r['n_pixels'] for r in primary):,}–"
                          f"{max(r['n_pixels'] for r in primary):,} jointly qualified pairs/band: bias (EnMAP minus S2) "
                          f"{min(r['bias_enmap_minus_actual_s2'] for r in primary):+.4f} to "
                          f"{max(r['bias_enmap_minus_actual_s2'] for r in primary):+.4f}; RMSD "
                          f"{min(r['rmsd'] for r in primary):.4f}–{max(r['rmsd'] for r in primary):.4f}. "
                          'One-day, atmosphere, view/sun, PSF and resampling differences prevent calibration-validation claims.')
    figure_note = ''
    for name in ('spectra_quality.png', 'transfer_summary.png'):
        if (root / cfg['figure_directory'] / name).exists():
            figure_note += f'\n![{name}]({cfg["figure_directory"]}/{name})\n'
    text = f'''# Real EnMAP first look — Berlin, 2025-05-11

**16-year calendar gap:** 2025 real EnMAP versus the 2009 simulated scene/reference.
These are historical discrepancies, combining model/sensor, geometry, season and
genuine land-cover change; they do not establish same-date transfer accuracy.

XML confirms centre acquisition **2025-05-11 10:42:03.869622 UTC**, processing-chain
version **{meta['processing_chain_revision']}**, 224 bands (418.416–2445.3 nm),
viewing zenith **{angle['viewingZenithAngle']['center']:.3f} degrees**,
viewing azimuth **{angle['viewingAzimuthAngle']['center']:.3f} degrees**,
solar elevation/azimuth **{angle['sunElevationAngle']['center']:.3f}/{angle['sunAzimuthAngle']['center']:.3f} degrees**.
Reflectance is **DN x 0.0001 + 0**, applied exactly once; raw nodata is **-32768**.
Land mode, bilinear, terrain correction, no cirrus/haze removal or band interpolation.

{quality_note}

{spectrum_note}

{metric_note}

Real spectra are bilinearly aligned to the exact tutorial grid and linearly
interpolated to its frozen 177-band axis; no extrapolation or gaps >35 nm.
At a gap boundary, nearest-band mapping is allowed only within half its FWHM and
5 nm (1783 maps to 1780.22 nm); the gap itself is never bridged.
This approximates differing spectral responses. No reference-fitted correction,
normalization, tuning or sum-to-one is applied. All 12 rebuilt EnMAP estimators
exactly reproduce saved Stage A predictions. The extra S2-like branch is a
sensitivity: one estimator differs at 2/71,158 replay pixels (maximum 0.009905).
Ensemble SD is across 3 runs, not calibrated uncertainty; nominal georeferencing
may contribute to historical differences.

{optional_note}

Raw files, model objects, detailed tables and fraction/uncertainty GeoTIFFs remain
local under ignored `data/real_enmap/20250511/`. Only code/configuration, this report
and two figures are published. Credentials are local and never logged. Deadline:
**2026-11-04 15:30 Europe/Berlin**. Re-run with `python real_enmap.py --stage metadata`,
`--stage download`, `--stage quality`, `--stage spectra`, `--stage transfer`.
Optional order: `--stage optional-prepare`, `optional-audit`, `optional-coreg`, `optional-compare`.

[DLR scene metadata](https://geoservice.dlr.de/eoc/ogc/stac/v1/collections/ENMAP_HSI_L2A/items/{cfg['scene_id']}) ·
[Quality codes: product specification section 4.5](https://www.enmap.org/data/doc/EN-PCV-ICD-2009-2_HSI_Product_Specification_Level1_Level2.pdf).
{figure_note}'''
    (root / cfg['report']).write_text(text, encoding='utf-8')
