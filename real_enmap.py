"""Real-data first-look stages; never tune the RF on the 2009 reference.

Raw data, local models and numerical diagnostics stay under ignored data/.
The same Stage A training seeds/settings are replayed because Stage A saved
ensemble predictions rather than serialized estimator objects.
"""
from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import joblib
import numpy as np
import yaml
from sklearn.ensemble import RandomForestRegressor

from mixing import mixtures
from real_enmap_io import Deadline, Downloader, StopWork, file_hash, utc_now, write_json


def parse_metadata(path):
    root = ET.parse(path).getroot()
    for element in root.iter():
        element.tag = element.tag.split('}')[-1]
    records = []
    for band in root.findall('specific/bandCharacterisation/bandID'):
        records.append({'band_1_based': int(band.attrib['number']),
                        'wavelength_nm': float(band.findtext('wavelengthCenterOfBand')),
                        'fwhm_nm': float(band.findtext('FWHMOfBand')),
                        'gain': float(band.findtext('GainOfBand')),
                        'offset': float(band.findtext('OffsetOfBand'))})
    records.sort(key=lambda band: band['band_1_based'])
    if not records or [band['band_1_based'] for band in records] != list(range(1, len(records) + 1)):
        raise StopWork('XML spectral band numbering is incomplete; no scaling is assumed.')
    if not all(np.isfinite([band['wavelength_nm'], band['fwhm_nm'], band['gain'], band['offset']]).all()
               and band['gain'] > 0 for band in records):
        raise StopWork('XML contains invalid spectral scaling or wavelength values.')
    angles = {}
    for name in ('viewingZenithAngle', 'viewingAzimuthAngle', 'acrossOffNadirAngle',
                 'alongOffNadirAngle', 'sunElevationAngle', 'sunAzimuthAngle'):
        parent = root.find('specific/' + name)
        angles[name] = {child.tag: float(child.text) for child in parent
                        if child.text and child.text.strip() != 'NA'} if parent is not None else None
    times = {point.findtext('frame'): point.findtext('utcTime')
             for point in root.findall('base/spatialCoverage/boundingPolygon/point')}
    revision = root.findtext('base/revision')
    background = root.findtext('specific/backgroundValue')
    if not revision or background is None:
        raise StopWork('XML is missing processing-chain revision or background value.')
    return {'source': 'downloaded METADATA.XML, not STAC substitution',
            'metadata_sha256': file_hash(path), 'processing_chain_revision': revision,
            'archived_L0_processor_version': root.findtext('base/archivedVersion'),
            'acquisition_times_utc_by_corner': times,
            'datatake_start_utc': root.findtext('specific/datatakeStart'),
            'datatake_stop_utc': root.findtext('specific/datatakeStop'),
            'angles_deg_by_corner': angles, 'background_raw_dn': float(background),
            'processing_parameters': {name: root.findtext('processing/' + name) for name in
                ('correctionType', 'cirrusHazeRemoval', 'bandInterpolation', 'imageResampling', 'terrainCorrection')},
            'physical_conversion': 'reflectance = raw_DN * per_band_gain + per_band_offset; applied exactly once',
            'band_count': len(records), 'bands': records}


def prepare_models(cfg, root, deadline):
    stage_a = yaml.safe_load((root / cfg['stage_a_config']).read_text(encoding='utf-8'))
    raw = root / cfg['raw_directory']
    raw.mkdir(parents=True, exist_ok=True)
    with np.load(root / stage_a['results'] / 'cache/prepared.npz') as prepared:
        library = prepared['library']
        labels = prepared['labels']
        wavelengths = prepared['wavelengths']
        weights = prepared['weights']
        valid = prepared['valid_image']
        pixels = prepared['image'][:, valid].T
    sensors = {'EnMAP': pixels, 'S2-like': pixels @ weights}
    model_path = raw / 'stage_a_models.joblib'
    models, audit = {}, []
    with np.load(root / stage_a['results'] / 'cache/ensemble_predictions.npz') as saved:
        for run in range(stage_a['mixing']['ensemble_size']):
            for target, name in enumerate(stage_a['classes']):
                deadline.check(10)
                seed = stage_a['seed'] + run * 100 + target
                x, y = mixtures(library, labels, target, stage_a['mixing'], seed)
                for sensor, simulation in sensors.items():
                    deadline.check(10)
                    model = RandomForestRegressor(**stage_a['forest'], random_state=seed)
                    model.fit(x if sensor == 'EnMAP' else x @ weights, y)
                    prediction = model.predict(simulation).astype('float32')
                    original = saved[sensor][run, target]
                    difference = float(np.max(np.abs(prediction - original)))
                    if sensor == 'EnMAP' and difference > 1e-7:
                        raise StopWork('Replayed model differs from saved Stage A predictions; transfer is not allowed.')
                    models[(sensor, run, target)] = model
                    audit.append({'sensor': sensor, 'member': run + 1, 'class': name,
                                  'seed': seed, 'max_abs_difference_from_saved_stage_a': difference,
                                  'pixels_changed_above_1e_7': int((np.abs(prediction - original) > 1e-7).sum()),
                                  'mean_abs_difference_from_saved_stage_a': float(np.abs(prediction - original).mean()),
                                  'exact_to_1e_7': difference <= 1e-7})
                    print(f'Stage A estimator replay verified: {sensor}, member {run + 1}, {name}', flush=True)
    joblib.dump({'models': models, 'wavelengths_nm': wavelengths, 's2_weights': weights,
                 'stage_a_config': stage_a}, model_path)
    result = {'time_utc': utc_now(), 'models_replayed': len(models),
              'maximum_absolute_prediction_difference': max(row['max_abs_difference_from_saved_stage_a'] for row in audit),
              'enmap_maximum_absolute_prediction_difference': max(row['max_abs_difference_from_saved_stage_a'] for row in audit if row['sensor'] == 'EnMAP'),
              's2_like_nonidentical_replay_is_sensitivity_only': any(not row['exact_to_1e_7'] for row in audit if row['sensor'] == 'S2-like'),
              'training_source': 'unchanged 75-spectrum tutorial library and synthetic mixtures',
              'input_units': 'Stage A raw reflectance scale (physical reflectance x 10000)',
              'no_reference_tuning': True, 'real_image_inference_completed': False,
              'audit': audit}
    write_json(raw / 'model_replay_audit.json', result)
    print('RF replay audit saved. EnMAP replay passed the strict baseline check; S2-like differences are explicitly recorded.', flush=True)
    return result


def run_metadata(cfg, root, deadline):
    downloader = Downloader(cfg, root, deadline)
    item = downloader.item()
    path = downloader.asset(item, 'metadata', 'metadata.xml')
    parsed = parse_metadata(path)
    write_json(root / cfg['raw_directory'] / 'metadata_summary.json', parsed)
    print(json.dumps({key: value for key, value in parsed.items() if key != 'bands'}, indent=2), flush=True)
    return parsed


def download_inputs(cfg, root, deadline):
    downloader = Downloader(cfg, root, deadline)
    item = downloader.item()
    keys = ['quality_classes', 'quality_cloud', 'quality_cloud_shadow', 'quality_haze',
            'quality_cirrus', 'quality_snow', 'quality_testflags', 'defective_pixel_mask', 'image']
    for key in keys:
        downloader.asset(item, key, key + '.tif')
    downloader.session.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='real_enmap_config.yaml')
    parser.add_argument('--stage', choices=['metadata', 'download', 'prepare-models', 'quality', 'spectra', 'transfer', 'optional-prepare', 'optional-audit', 'optional-coreg', 'optional-compare'], required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    cfg = yaml.safe_load((root / args.config).read_text(encoding='utf-8'))
    deadline = Deadline(cfg['deadline'], root, cfg['report'])
    try:
        deadline.arm()
        if args.stage == 'metadata':
            run_metadata(cfg, root, deadline)
        elif args.stage == 'download':
            download_inputs(cfg, root, deadline)
        elif args.stage in ['quality', 'spectra', 'transfer']:
            from real_enmap_processing import quality_stage, spectra_stage, transfer_stage
            {'quality': quality_stage, 'spectra': spectra_stage, 'transfer': transfer_stage}[args.stage](cfg, root, deadline)
        elif args.stage in ['optional-prepare', 'optional-audit', 'optional-coreg', 'optional-compare']:
            from real_enmap_optional import prepare_s2, audit_s2_radiometry, run_coreg, compare_reflectance
            {'optional-prepare': prepare_s2, 'optional-audit': audit_s2_radiometry, 'optional-coreg': run_coreg,
             'optional-compare':compare_reflectance}[args.stage](cfg, root, deadline)
        else:
            prepare_models(cfg, root, deadline)
    except StopWork as exc:
        write_json(root / cfg['raw_directory'] / 'last_stop.json',
                   {'time_utc': utc_now(), 'stage': args.stage, 'reason': str(exc)})
        print('STOP: ' + str(exc), flush=True)
        raise SystemExit(3)
    finally:
        deadline.close()


if __name__ == '__main__':
    main()
