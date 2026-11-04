"""Independent external-endmember preflight and gated spectral-temporal study."""
from __future__ import annotations

import json
import logging
import math
import time
import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from affine import Affine
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform_bounds

from data_io import write_raster

LOG = logging.getLogger(__name__)
WORLD_COVER_CODES = {50: 1, 60: 3, 80: 4, 10: 2, 20: 2, 30: 2, 40: 2, 90: 2, 95: 2, 100: 2}


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)), encoding='utf-8')


def stac_search(endpoint, collection, bbox, interval, cache):
    request_key={'endpoint':endpoint,'collection':collection,'bbox':list(bbox),'datetime':interval}
    if cache.exists():
        stored=json.loads(cache.read_text(encoding='utf-8'))
        if stored.get('_request')==request_key:
            return stored['features']
    # STAC is an explicit public API, not HTML scraping.
    url = endpoint.rstrip('/')+'/search'
    body = {'collections': [collection], 'bbox': list(bbox), 'datetime': interval, 'limit': 100}
    features = []
    while url:
        response = requests.post(url, json=body, timeout=60) if body is not None else requests.get(url, timeout=60)
        response.raise_for_status()
        result = response.json()
        features.extend(result.get('features', []))
        link = next((link for link in result.get('links', []) if link.get('rel') == 'next'), None)
        if link is None:
            break
        url = link['href']
        if link.get('method', 'GET').upper() == 'POST':
            body = ({**body, **link.get('body', {})} if link.get('merge') and body else link.get('body'))
        else:
            body = None
    dump(cache, {'type': 'FeatureCollection', 'features': features,'_request':request_key})
    return features


def stop(out, reason, evidence):
    dump(out/'gate_status.json', {'status': 'stopped', 'reason': reason, 'evidence': evidence})
    (out/'STATUS.md').write_text(
        '# Stage B stopped at a fixed gate\n\n'+reason+'\n\n'+
        'Stage A is complete and unchanged. Temporal/fusion model results and a value-of-acquisition map are not claimed.\n\n'+
        'The 2009 simulated EnMAP, 2017/2018 S2 candidates, and 2021 WorldCover are not contemporaneous.\n\n'+
        'Evidence is in `gate_status.json`, the WorldCover audits and cached STAC metadata.\n', encoding='utf-8')
    LOG.warning('Stage B stopped: %s', reason)
    return {'status': 'stopped', 'reason': reason}


def worldcover_preflight(cfg, root, out, profile, reference, valid):
    ext = cfg['extension']
    pad = int(math.ceil(ext['buffer_outer_m']/30))
    h, w = profile['height'], profile['width']
    shape = (h+2*pad, w+2*pad)
    tx = profile['transform'] * Affine.translation(-pad, -pad)
    bounds = rasterio.transform.array_bounds(*shape, tx)
    bbox = transform_bounds(profile['crs'], 'EPSG:4326', *bounds)
    cache = out/'cube'
    cache.mkdir(exist_ok=True)
    year = ext['worldcover_year']
    items = stac_search('https://planetarycomputer.microsoft.com/api/stac/v1', 'esa-worldcover', bbox,
                        f'{year}-01-01T00:00:00Z/{year}-12-31T23:59:59Z', cache/'worldcover_stac.json')
    if not items:
        raise RuntimeError('WorldCover STAC returned no coverage for the external buffer.')
    grid_key=hashlib.sha256(json.dumps({'year':year,'shape':list(shape),'transform':list(tx)[:6],
                                      'crs':str(profile['crs'])},sort_keys=True).encode()).hexdigest()[:12]
    mapped_path = cache/f'worldcover_{year}_{grid_key}_10m.npy'
    sources = []
    if mapped_path.exists():
        fine = np.load(mapped_path)
    else:
        fine = np.zeros((shape[0]*3, shape[1]*3), dtype='uint8')
        for item in items:
            # The identical ESA product is public on AWS; avoid expiring signed
            # Azure URLs in reproducibility artifacts. Keep the STAC origin too.
            product = item['id']+'_Map.tif'
            version=re.search(r'_v(\d+)_',item['id']).group(1)
            url = f'https://esa-worldcover.s3.eu-central-1.amazonaws.com/v{version}/{year}/map/{product}'
            sources.append({'id': item['id'], 'stac_asset': item['assets']['map']['href'], 'public_esa_asset': url})
            LOG.info('Read independent WorldCover tile %s', item['id'])
            with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', GDAL_HTTP_TIMEOUT=60, GDAL_HTTP_MAX_RETRY=2):
                with rasterio.open(url) as source:
                    with WarpedVRT(source, crs=profile['crs'], transform=tx*Affine.scale(1/3, 1/3),
                                   width=shape[1]*3, height=shape[0]*3, resampling=Resampling.nearest,
                                   nodata=0, dtype='uint8') as vrt:
                        values = vrt.read(1)
                        np.copyto(fine, values, where=(fine == 0) & (values > 0))
        np.save(mapped_path, fine)
    if not sources:
        for item in items:
            version=re.search(r'_v(\d+)_',item['id']).group(1)
            sources.append({'id': item['id'], 'stac_asset': item['assets']['map']['href'],
                            'public_esa_asset': f'https://esa-worldcover.s3.eu-central-1.amazonaws.com/v{version}/{year}/map/{item["id"]}_Map.tif'})
    dump(out/'worldcover_sources.json', {'year': year, 'licence': 'CC-BY-4.0', 'sources': sources,
                                        'crosswalk': WORLD_COVER_CODES,
                                        'warning': 'WorldCover labels are not physical endmember truth and are 12 years later than the 2009 scene.'})
    grouped = np.zeros_like(fine)
    for source_code, mapped_code in WORLD_COVER_CODES.items():
        grouped[fine == source_code] = mapped_code
    fractions = np.asarray([(grouped == code).reshape(shape[0], 3, shape[1], 3).mean(axis=(1, 3))
                            for code in range(1, 5)], dtype='float32')
    crop = fractions[:, pad:pad+h, pad:pad+w]
    dominant = crop.argmax(axis=0)
    known = crop.sum(axis=0) > 0
    stable = known & (dominant == reference.argmax(axis=0)) & valid
    agreement = float(stable.sum()/valid.sum())
    write_raster(out/'worldcover_fractions_2021.tif', crop, profile, cfg['classes'])
    write_raster(out/'worldcover_reference_consistency_mask.tif', np.where(valid, stable, np.nan)[None], profile, ['worldcover_reference_consistency'])
    confusion = pd.crosstab(pd.Series(reference[:, valid].argmax(axis=0), name='reference'),
                            pd.Series(np.where(known, dominant, -1)[valid], name='worldcover'))
    confusion.to_csv(out/'worldcover_reference_confusion.csv')
    yy, xx = np.indices(shape)
    xs = tx.c + (xx+.5)*30
    ys = tx.f - (yy+.5)*30
    left, bottom, right, top = rasterio.transform.array_bounds(h, w, profile['transform'])
    dx = np.maximum(np.maximum(left-xs, xs-right), 0)
    dy = np.maximum(np.maximum(bottom-ys, ys-top), 0)
    distance = np.hypot(dx, dy)
    outside = (distance > ext['buffer_inner_m']) & (distance <= ext['buffer_outer_m'])
    rng = np.random.default_rng(cfg['seed'])
    counts, candidates = [], []
    for target, name in enumerate(cfg['classes']):
        r, c = np.where(outside & (fractions[target] >= ext['minimum_endmember_purity']))
        n = len(r)
        take = rng.choice(n, size=min(n, ext['maximum_endmembers_per_class']), replace=False)
        counts.append({'class': name, 'pure_external_pixels_available': n, 'sampled_candidates': len(take), 'worldcover_year': year})
        for index in take:
            candidates.append({'class': name, 'class_index': target, 'row': int(r[index]), 'col': int(c[index]),
                               'x': float(xs[r[index], c[index]]), 'y': float(ys[r[index], c[index]]),
                               'worldcover_purity': float(fractions[target, r[index], c[index]]), 'worldcover_year': year})
    pd.DataFrame(counts).to_csv(out/'endmember_support_preflight.csv', index=False)
    pd.DataFrame(candidates).to_csv(out/'external_endmember_candidates.csv', index=False)
    evidence = {'worldcover_year': year, 'reference_pixels': int(valid.sum()),
                'worldcover_reference_agreement': agreement, 'consistent_pixels': int(stable.sum()),
                'inconsistent_pixels': int(valid.sum()-stable.sum()), 'endmember_counts': counts,
                'buffer_inner_m': ext['buffer_inner_m'], 'buffer_outer_m': ext['buffer_outer_m']}
    dump(out/'worldcover_preflight.json', evidence)
    return evidence, stable, pd.DataFrame(candidates), tx, shape


def run_extension(cfg, root, results, profile, reference, valid, preflight_only=False):
    out = results/'stage_b'
    out.mkdir(exist_ok=True)
    history_path=out/'runtime_history.json'
    history=json.loads(history_path.read_text(encoding='utf-8')) if history_path.exists() else []
    if not history and (out/'runtime.json').exists():
        earlier=json.loads((out/'runtime.json').read_text(encoding='utf-8'))
        prior_status=json.loads((out/'gate_status.json').read_text(encoding='utf-8')).get('status','unknown') if (out/'gate_status.json').exists() else 'unknown'
        history.append({'seconds':earlier['seconds'],'status':prior_status,'note':'Earlier download/preprocessing attempt, retained before cached retry.'})
    import yaml
    (out/'config_used.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    dump(out/'gate_status.json', {'status':'running','phase':'worldcover_preflight'})
    start = time.perf_counter()
    try:
        evidence, stable, candidates, tx, shape = worldcover_preflight(cfg, root, out, profile, reference, valid)
        missing = [row['class'] for row in evidence['endmember_counts']
                   if row['pure_external_pixels_available'] < cfg['extension']['minimum_endmembers_per_class']]
        if missing:
            return stop(out, f'External WorldCover endmember support is below 30 for: {missing}. Purity and class definitions were not relaxed.', evidence)
        if evidence['worldcover_reference_agreement'] < cfg['extension']['minimum_worldcover_reference_agreement']:
            return stop(out, f"WorldCover/reference dominant-class agreement is {evidence['worldcover_reference_agreement']:.1%}, below the fixed 70% gate.", evidence)
        if preflight_only:
            LOG.info('WorldCover support and agreement gates passed: %s', evidence)
            return {'status': 'worldcover_preflight_passed', 'evidence': evidence}
        dump(out/'gate_status.json', {'status':'running','phase':'temporal_observations_and_features','evidence':evidence})
        from temporal import temporal_analysis
        return temporal_analysis(cfg, root, out, profile, reference, valid, evidence, stable, candidates, tx, shape)
    except Exception as exc:
        dump(out/'gate_status.json', {'status': 'technical_failure', 'error_type': type(exc).__name__, 'reason': str(exc)})
        (out/'STATUS.md').write_text('# Stage B technical failure\n\n'+str(exc)+'\n\nStage A remains complete. This is an access/software failure, not a scientific gate result.\n', encoding='utf-8')
        LOG.exception('Stage B technical failure')
        raise
    finally:
        seconds=time.perf_counter()-start
        status=json.loads((out/'gate_status.json').read_text(encoding='utf-8')).get('status','unknown')
        history.append({'seconds':seconds,'status':status})
        dump(history_path,history)
        dump(out/'runtime.json', {'seconds':seconds,'total_recorded_wall_seconds':sum(row['seconds'] for row in history),
                                 'note':'Separate attempts include initial remote downloads and subsequent cached processing.'})
