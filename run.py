"""Reproduce all configured stages: python run.py --config config.yaml."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import logging
import os
import platform
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for key, value in {'MPLBACKEND': 'Agg', 'MPLCONFIGDIR': str(ROOT/'.tmp/mpl'),
                   'TEMP': str(ROOT/'.tmp'), 'TMP': str(ROOT/'.tmp'),
                   'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}.items():
    os.environ.setdefault(key, value)
(ROOT/'.tmp').mkdir(exist_ok=True)

import numpy as np
import pandas as pd
import yaml

from data_io import load_inputs, reference_fractions, sensor_weights, write_raster, sha256
from evaluation import differences, summarize
from mixing import train_sensor_pair


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--stage', choices=['prepare', 'a', 'a-additions', 'b', 'all'], default='all')
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    cfg = yaml.safe_load(config_path.read_text(encoding='utf-8'))
    root = config_path.parent
    out = root / cfg['results']
    out.mkdir(exist_ok=True)
    (out/'cache').mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out/'run.log', encoding='utf-8')])
    logging.getLogger('rasterio.session').setLevel(logging.WARNING)
    log = logging.getLogger(__name__)
    start = time.perf_counter()
    snapshot_name = {'b':'stage_b_config_used.yaml', 'a-additions':'stage_a_additions_config_used.yaml'}.get(args.stage,'config_used.yaml')
    (out/snapshot_name).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    environment = {'python': sys.version, 'platform': platform.platform(),
                   'packages': {d.metadata['Name']: d.version for d in importlib.metadata.distributions()}}
    (out/'environment.json').write_text(json.dumps(environment, indent=2), encoding='utf-8')
    if args.stage == 'a-additions':
        from stage_a_additions import run_additions
        run_additions(cfg, root, out)
        return
    if args.stage == 'b' and (out/'cache/prepared.npz').exists():
        previous = yaml.safe_load((out/'config_used.yaml').read_text(encoding='utf-8'))
        if any(previous[key] != cfg[key] for key in ['data', 'classes', 'reference', 'spectral']):
            raise ValueError('Stage A input/reference settings changed; re-run prepare before Stage B.')
        manifest = json.loads((out/'input_manifest.json').read_text(encoding='utf-8'))
        for name, detail in manifest['files'].items():
            if sha256(root/cfg['data']['directory']/name) != detail['sha256']:
                raise ValueError(f'Input changed since preparation: {name}')
        with np.load(out/'cache/prepared.npz') as prepared:
            reference, valid = prepared['reference'], prepared['valid']
        import rasterio
        with rasterio.open(root/cfg['data']['directory']/'enmap_berlin.tif') as src:
            profile = src.profile.copy()
        from stage_b import run_extension
        run_extension(cfg, root, out, profile, reference, valid)
        return
    image, library, labels, wavelengths, valid_image, profile, manifest = load_inputs(cfg, root, out)
    reference, valid, polygons = reference_fractions(cfg, root, profile, valid_image, out)
    weights, srf_info = sensor_weights(cfg, root, wavelengths, out)
    np.savez_compressed(out/'cache/prepared.npz', image=image, library=library, labels=labels,
                        wavelengths=wavelengths, valid_image=valid_image, reference=reference, valid=valid, weights=weights)
    if args.stage == 'prepare':
        log.info('Preparation complete: %.2f seconds', time.perf_counter()-start)
        return
    if args.stage in ['a', 'all']:
        predictions, importance, timings = train_sensor_pair(cfg, image, library, labels, valid_image, weights, out)
        summaries, block_tables, context_tables, difference_tables, interval_tables = [], [], [], [], []
        baseline_means = baseline_stds = None
        h, w = valid.shape
        for size in sorted(set([cfg['mixing']['ensemble_size'], cfg['mixing']['sensitivity_ensemble_size']])):
            means, stds = {}, {}
            for sensor, ensemble in predictions.items():
                means[sensor] = np.full((len(cfg['classes']), h, w), np.nan, dtype='float32')
                stds[sensor] = means[sensor].copy()
                means[sensor][:, valid_image] = ensemble[:size].mean(axis=0)
                stds[sensor][:, valid_image] = ensemble[:size].std(axis=0, ddof=0)
                stem = sensor.lower().replace('-', '_')
                write_raster(out/f'{stem}_fractions_n{size}.tif', means[sensor], profile, cfg['classes'])
                write_raster(out/f'{stem}_uncertainty_n{size}.tif', stds[sensor], profile, cfg['classes'])
            summary, blocks, contexts = summarize(cfg, reference, valid, means, stds, size)
            diff, intervals = differences(cfg, reference, valid, means, blocks, size)
            summaries.append(summary); block_tables.append(blocks); context_tables.append(contexts)
            difference_tables.append(diff); interval_tables.append(intervals)
            if size == cfg['mixing']['ensemble_size']:
                baseline_means, baseline_stds = means, stds
        for name, tables in [('summary', summaries), ('blocks', block_tables), ('contexts', context_tables),
                             ('block_differences', difference_tables), ('paired_block_bootstrap', interval_tables)]:
            pd.concat(tables, ignore_index=True).to_csv(out/f'{name}.csv', index=False)
        pd.DataFrame(timings).to_csv(out/'model_runtime.csv', index=False)
        for sensor, values in importance.items():
            axis = wavelengths if sensor == 'EnMAP' else srf_info.effective_wavelength_nm.to_numpy()
            table = pd.DataFrame(values[:cfg['mixing']['ensemble_size']].mean(axis=0).T, columns=cfg['classes'])
            table.insert(0, 'wavelength_nm', axis)
            table.to_csv(out/f'{sensor.lower().replace("-", "_")}_feature_importance.csv', index=False)
        from plotting import figures
        figures(cfg, out, reference, valid, profile, baseline_means, baseline_stds, wavelengths,
                importance['EnMAP'][:cfg['mixing']['ensemble_size']].mean(axis=0))
        from report import write_report
        elapsed = time.perf_counter()-start
        (out/'runtime.json').write_text(json.dumps({'stage_a_seconds': elapsed, 'model_seconds_sum': sum(t['seconds'] for t in timings),
                                                   'note': 'Wall time excludes initial environment installation and input downloads.'}, indent=2), encoding='utf-8')
        write_report(cfg, root, out, manifest, srf_info, elapsed)
        log.info('Stage A finished in %.1f seconds', elapsed)
        if cfg.get('stage_a_additions',{}).get('enabled',False):
            from stage_a_additions import run_additions
            run_additions(cfg, root, out)
    if args.stage == 'b' or (args.stage == 'all' and cfg['extension']['enabled']):
        from stage_b import run_extension
        run_extension(cfg, root, out, profile, reference, valid)


if __name__ == '__main__':
    main()
