"""Optional actual-S2 comparison, after the four required first-look stages.

AROSICS is run in an isolated conda-forge environment with a short Windows path.
No 2009 reference labels or unmixing errors are used for matching or calibration.
"""
from __future__ import annotations

import json
from pathlib import Path

BANDS = ['blue', 'green', 'red', 'rededge1', 'rededge2', 'rededge3', 'nir', 'nir08', 'swir16', 'swir22']
S2_BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']


def asset_on_grid(asset, target):
    import numpy as np
    import rasterio
    from rasterio.warp import reproject, Resampling, transform_bounds
    from real_enmap_processing import crop_window
    href = asset['href']
    if href.startswith('s3://'):
        bucket, _, key = href[5:].partition('/')
        href = 'https://' + bucket + '.s3.eu-central-1.amazonaws.com/' + key
    bounds = rasterio.transform.array_bounds(target['height'], target['width'], target['transform'])
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', GDAL_HTTP_TIMEOUT='30',
                      GDAL_HTTP_MAX_RETRY='1', AWS_NO_SIGN_REQUEST='YES'):
        with rasterio.open(href) as src:
            source_bounds = transform_bounds(target['crs'], src.crs, *bounds, densify_pts=21)
            window = crop_window(src, source_bounds)
            data = src.read(1, window=window).astype('float32')
            source_transform, source_crs = src.window_transform(window), src.crs
    result = np.full((target['height'], target['width']), np.nan, dtype='float32')
    reproject(np.where(data != 0, data, np.nan), result, src_transform=source_transform, src_crs=source_crs,
              dst_transform=target['transform'], dst_crs=target['crs'], src_nodata=np.nan, dst_nodata=np.nan,
              resampling=Resampling.average, num_threads=2)
    return result, data, source_transform, source_crs


def prepare_s2(cfg, root, deadline):
    import numpy as np
    import rasterio
    from affine import Affine
    from pyproj import Transformer
    from shapely.geometry import shape, box
    from shapely import transform
    from real_enmap_io import StopWork, write_json, utc_now
    from real_enmap_processing import warp, load_context
    from data_io import write_raster
    raw, meta, prepared, profile, bounds = load_context(cfg, root)
    if not (raw / 'transfer_summary.json').exists():
        raise StopWork('Optional work requires completed required first-look stages.')
    deadline.check(60)
    search = json.loads((raw / 's2_search.json').read_text(encoding='utf-8'))
    item = next(f for f in search['features'] if f['id'] == cfg['optional_cross_mission']['selected_s2_item'])
    to_projected = Transformer.from_crs(4326, profile['crs'], always_xy=True).transform
    coverage = transform(shape(item['geometry']), to_projected, interleaved=False).intersection(box(*bounds)).area / box(*bounds).area
    if coverage < .99:
        raise StopWork('Selected S2 item covers less than 99% of the tutorial; do not silently compare another area.')
    scl_grid, scl, scl_transform, scl_crs = asset_on_grid(item['assets']['scl'], profile)
    scl_profile = profile.copy(); scl_profile.update(transform=scl_transform, crs=scl_crs, height=scl.shape[0], width=scl.shape[1])
    scl_clear = np.isin(scl, [4, 5, 6])
    clear = warp(scl_clear, scl_profile, profile, rasterio.warp.Resampling.average) >= 1 - 1e-6
    values, raw_grids, calibration = [], [], []
    for band in BANDS:
        deadline.check(30)
        scaled_dn, _, _, _ = asset_on_grid(item['assets'][band], profile)
        calibration_meta = item['assets'][band]['raster:bands'][0]
        scale, offset = calibration_meta['scale'], calibration_meta.get('offset', 0)
        raw_grids.append(scaled_dn)
        values.append(scaled_dn * scale + offset)
        calibration.append({'band': band, 'scale': scale, 'declared_offset': offset,
                            'effective_offset': offset, 'provider_boa_offset_applied': item['properties'].get('earthsearch:boa_offset_applied')})
        print('Actual S2 band prepared: ' + band, flush=True)
    actual = np.asarray(values, dtype='float32')
    raw_grids = np.asarray(raw_grids, dtype='float32')
    negative = {band: float((actual[i][clear] < 0).mean()) for i, band in enumerate(BANDS)}
    np.savez_compressed(raw / 's2_actual_prepared.npz', reflectance=actual, raw_dn=raw_grids, clear=clear)
    with np.load(raw / 'quality_prepared.npz') as cache:
        native = cache['source_diagnostic_reflectance']; axis = cache['source_axis_nm']
        native_transform = Affine(*cache['native_transform'])
    near = int(np.argmin(abs(axis - 842)))
    native_profile = profile.copy(); native_profile.update(transform=native_transform, height=native.shape[1], width=native.shape[2])
    enmap_nir = warp(native[near], native_profile, profile, fill=np.nan)
    good = clear & np.isfinite(actual[6]) & np.isfinite(enmap_nir)
    write_raster(raw / 'coreg_s2_nir.tif', np.where(clear & np.isfinite(actual[6]), actual[6], np.nan)[None], profile, ['actual_S2C_B8'])
    write_raster(raw / 'coreg_enmap_nir.tif', enmap_nir[None], profile, ['EnMAP_nearest_842nm'])
    points = []
    ws = cfg['optional_cross_mission']['coreg_window_pixels']
    for row in [60, 100, 140, 180, 220, 260, 300, 350]:
        for column in [45, 75, 110, 145, 175]:
            if good[row, column]:
                fraction = float(good[row-ws//2:row+ws//2, column-ws//2:column+ws//2].mean())
                points.append((fraction, row, column))
    points.sort(reverse=True)
    positions = []
    for fraction, row, column in points:
        if fraction < .8:
            continue
        if any(abs(row - q['row']) < 64 and abs(column - q['column']) < 64 for q in positions):
            continue
        x, y = profile['transform'] * (column + .5, row + .5)
        positions.append({'row':row,'column':column,'x':x,'y':y,'joint_good_fraction':fraction})
        if len(positions) == 4:
            break
    write_json(raw / 'optional_geometry_inputs.json', {'raw_directory': str(raw), 'points': positions,
               'window_size': ws, 'max_shift_pixels': cfg['optional_cross_mission']['coreg_max_shift_pixels'],
               'minimum_reliability_percent': cfg['optional_cross_mission']['coreg_min_reliability_percent']})
    result = {'time_utc': utc_now(), 'item_id': item['id'], 'acquisition_datetime': item['properties']['datetime'],
              'day_gap_from_enmap': 1, 'scene_cloud_percent': item['properties']['eo:cloud_cover'],
              'footprint_coverage': float(coverage), 'platform': item['properties']['platform'],
              'conservative_scl_clear_pixels': int(clear.sum()), 'joint_nir_good_pixels': int(good.sum()),
              'calibration': calibration, 'negative_reflectance_share_by_band': negative,
              'radiometric_rule': 'Declared per-asset raster:bands scale/offset, exactly once, following Earth Search producer documentation.'}
    write_json(raw / 's2_actual_summary.json', result)
    print(json.dumps(result, indent=2), flush=True)
    return result


def coreg_worker(raw):
    import numpy as np
    import arosics
    from arosics import COREG
    raw = Path(raw)
    cfg = json.loads((raw / 'optional_geometry_inputs.json').read_text(encoding='utf-8'))
    rows = []
    for point in cfg['points']:
        try:
            model = COREG(str(raw/'coreg_s2_nir.tif'), str(raw/'coreg_enmap_nir.tif'),
                          wp=(point['x'], point['y']), ws=(cfg['window_size'], cfg['window_size']),
                          max_shift=cfg['max_shift_pixels'], CPUs=1, progress=False, q=True, nodata=(-9999,-9999))
            model.calculate_spatial_shifts()
            reliability = float(model.shift_reliability)
            row = dict(point, success=bool(model.success), reliability_percent=reliability,
                       x_correction_m=float(model.x_shift_map), y_correction_m=float(model.y_shift_map))
            row['accepted'] = row['success'] and reliability >= cfg['minimum_reliability_percent']
        except Exception as error:
            row = dict(point, success=False, accepted=False, error_type=type(error).__name__)
        rows.append(row)
    accepted = [row for row in rows if row['accepted']]
    result = {'arosics_version':arosics.__version__, 'reference':'actual S2C 2025-05-12 B8',
              'target':'real EnMAP 2025-05-11 nearest 842nm band', 'windows':rows,
              'input_mask_rule':cfg.get('mask_rule','unspecified'),
              'accepted_windows':len(accepted), 'shift_applied_to_required_model_results':False}
    if accepted:
        result['median_x_correction_m'] = float(np.median([row['x_correction_m'] for row in accepted]))
        result['median_y_correction_m'] = float(np.median([row['y_correction_m'] for row in accepted]))
        result['correction_direction'] = 'Add these map-coordinate corrections to EnMAP origin to align it to S2.'
    (raw/'arosics_result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


def audit_s2_radiometry(cfg, root, deadline):
    import numpy as np
    import rasterio
    import requests
    import xml.etree.ElementTree as ET
    from real_enmap_io import StopWork, write_json, utc_now
    from data_io import write_raster
    raw = root / cfg['raw_directory']
    search = json.loads((raw/'s2_search.json').read_text(encoding='utf-8'))
    item = next(f for f in search['features'] if f['id']==cfg['optional_cross_mission']['selected_s2_item'])
    with np.load(raw/'s2_actual_prepared.npz') as z:
        dn, clear = z['raw_dn'], z['clear']
    with rasterio.open(root/'data/enmap_berlin.tif') as src:
        profile = src.profile.copy()
    meta_url = item['assets']['product_metadata']['href'].replace('s3://sentinel-s2-l2a/','https://sentinel-s2-l2a.s3.eu-central-1.amazonaws.com/')
    session = requests.Session(); session.auth = lambda request: request
    response = session.get(meta_url,timeout=45)
    if response.status_code == 404:
        info_url = item['assets']['tileinfo_metadata']['href'].replace('s3://sentinel-s2-l2a/','https://sentinel-s2-l2a.s3.eu-central-1.amazonaws.com/')
        info_response = session.get(info_url,timeout=45); info_response.raise_for_status()
        product_path = info_response.json()['productPath']
        if not product_path.startswith('products/') or '..' in product_path.split('/'):
            raise StopWork('Unexpected original Sentinel product metadata path.')
        meta_url = 'https://sentinel-s2-l2a.s3.eu-central-1.amazonaws.com/' + product_path + '/metadata.xml'
        response = session.get(meta_url,timeout=45)
    response.raise_for_status()
    xml = ET.fromstring(response.content)
    offsets = [float(e.text) for e in xml.iter() if e.tag.split('}')[-1]=='BOA_ADD_OFFSET']
    quantification = [float(e.text) for e in xml.iter() if e.tag.split('}')[-1]=='BOA_QUANTIFICATION_VALUE']
    if not offsets or len(set(offsets)) != 1 or quantification != [10000.0]:
        raise StopWork('Original S2 product metadata has unexpected radiometric conventions; no offset is guessed.')
    comparisons = []
    for band in ['red','nir','swir22']:
        deadline.check(90)
        original,_,_,_ = asset_on_grid(item['assets'][band+'-jp2'],profile)
        cog = dn[BANDS.index(band)]
        valid = clear & np.isfinite(original) & np.isfinite(cog) & (original>1000) & (cog>0)
        delta = cog[valid]-original[valid]
        if len(delta)<100:
            raise StopWork('Too few independent JP2/COG audit pixels.')
        median = float(np.median(delta))
        residual95 = float(np.quantile(abs(delta-median),.95))
        comparisons.append({'band':band,'n_pixels':int(valid.sum()),'median_cog_minus_original_jp2_dn':median,
                            'p95_absolute_departure_from_median_dn':residual95})
        print('Independent original JP2 audit:',band,median,'DN',flush=True)
    delta = float(np.median([row['median_cog_minus_original_jp2_dn'] for row in comparisons]))
    if any(abs(row['median_cog_minus_original_jp2_dn']-delta)>2 or
           row['p95_absolute_departure_from_median_dn']>2 for row in comparisons):
        raise StopWork('COG/JP2 differences are not a consistent offset; no reflectance correction is fitted.')
    if min(abs(delta),abs(delta-offsets[0]))>2:
        raise StopWork('COG/JP2 shift does not match the original documented BOA offset.')
    physical_offset = offsets[0]/quantification[0] - delta/quantification[0]
    # Round the independently established exact integer-DN convention, not image errors.
    physical_offset = round(physical_offset,4)
    actual = dn * .0001 + physical_offset
    np.savez_compressed(raw/'s2_actual_prepared.npz',reflectance=actual,raw_dn=dn,clear=clear)
    summary = json.loads((raw/'s2_actual_summary.json').read_text(encoding='utf-8'))
    for row in summary['calibration']:
        row['effective_offset'] = physical_offset
        row['independent_original_product_audit'] = True
    summary['radiometric_rule'] = 'Original ESA JP2/COG DN audit plus original BOA_ADD_OFFSET/QUANTIFICATION; no EnMAP or reference-label fitting.'
    summary['negative_reflectance_share_by_band'] = {band:float((actual[i][clear]<0).mean()) for i,band in enumerate(BANDS)}
    write_json(raw/'s2_actual_summary.json',summary)
    good = clear & np.isfinite(actual[6])
    write_raster(raw/'coreg_s2_nir.tif',np.where(good,actual[6],np.nan)[None],profile,['actual_S2C_B8_original_product_audited'])
    result = {'time_utc':utc_now(),'item_id':item['id'],'original_boa_add_offset_dn':offsets[0],
              'original_quantification':quantification[0],'original_boa_offset_band_count':len(offsets),
              'comparisons':comparisons,'effective_cog_physical_offset':physical_offset,
              'no_enmap_reference_or_prediction_errors_used':True,
              'declared_asset_offset':item['assets']['red']['raster:bands'][0].get('offset'),
              'explanation':'COG DN differs from the producer JP2 by the documented offset; applying that offset twice is prevented.'}
    write_json(raw/'s2_original_radiometry_audit.json',result)
    print(json.dumps(result,indent=2),flush=True)
    return result


def refresh_geometry_inputs(cfg, root, deadline):
    """Build independent source masks; do not give both FFT inputs identical holes."""
    import numpy as np
    import rasterio
    from affine import Affine
    from data_io import write_raster
    from real_enmap_processing import load_context, warp
    from real_enmap_io import write_json, StopWork
    raw,meta,prepared,profile,bounds=load_context(cfg,root)
    with np.load(raw/'s2_actual_prepared.npz') as z:
        s2=z['reflectance'][6];s2_clear=z['clear']
    with np.load(raw/'quality_prepared.npz') as z:
        native=z['source_diagnostic_reflectance'];axis=z['source_axis_nm'];tr=Affine(*z['native_transform'])
    native_profile=profile.copy();native_profile.update(transform=tr,height=native.shape[1],width=native.shape[2])
    enmap=warp(native[int(np.argmin(abs(axis-842)))],native_profile,profile,fill=np.nan)
    good=s2_clear&np.isfinite(s2)&np.isfinite(enmap)
    write_raster(raw/'coreg_s2_nir.tif',np.where(s2_clear,s2,np.nan)[None],profile,['S2C_B8_independent_mask'])
    write_raster(raw/'coreg_enmap_nir.tif',enmap[None],profile,['EnMAP_842_independent_mask'])
    ws=cfg['optional_cross_mission']['coreg_window_pixels'];candidates=[]
    for row in [60,100,140,180,220,260,300,350]:
        for column in [45,75,110,145,175]:
            if good[row,column]:
                fraction=float(good[row-ws//2:row+ws//2,column-ws//2:column+ws//2].mean())
                candidates.append((fraction,row,column))
    candidates.sort(reverse=True);points=[]
    for fraction,row,column in candidates:
        if fraction<.8:continue
        if any(abs(row-p['row'])<64 and abs(column-p['column'])<64 for p in points):continue
        x,y=profile['transform']*(column+.5,row+.5)
        points.append({'row':row,'column':column,'x':x,'y':y,'joint_good_fraction':fraction})
        if len(points)==4:break
    if not points:raise StopWork('No >=80% joint-good 64-pixel matching windows; no shift is guessed.')
    initial=raw/'arosics_result.json'
    if initial.exists() and not (raw/'arosics_initial_joint_mask_result.json').exists():
        (raw/'arosics_initial_joint_mask_result.json').write_bytes(initial.read_bytes())
    write_json(raw/'optional_geometry_inputs.json',{'raw_directory':str(raw),'points':points,
               'window_size':ws,'max_shift_pixels':cfg['optional_cross_mission']['coreg_max_shift_pixels'],
               'minimum_reliability_percent':cfg['optional_cross_mission']['coreg_min_reliability_percent'],
               'mask_rule':'Independent per-source quality masks; matching centres selected only by joint QA coverage.'})
    print('Independent-mask matching windows:',points,flush=True)


def run_coreg(cfg, root, deadline):
    import os
    import subprocess
    from real_enmap_io import StopWork, file_hash, write_json
    deadline.check(600)
    raw = root / cfg['raw_directory']
    refresh_geometry_inputs(cfg,root,deadline)
    tools = raw / 'optional_tools'
    status = json.loads((tools / 'installation_status.json').read_text(encoding='utf-8'))
    if status['state'] != 'installed':
        raise StopWork('AROSICS environment was not installed; no substitute algorithm is labelled AROSICS.')
    environment = os.environ.copy()
    environment['MAMBA_ROOT_PREFIX'] = str(root.parent / '.ar_root')
    command = [str(tools/'micromamba.exe'), '--no-rc', 'run', '--prefix', str(root.parent/'.ar_env'),
               'python', str(root/'real_enmap_optional.py'), '--coreg-worker-raw', str(raw)]
    with (raw/'arosics_run.log').open('w',encoding='utf-8') as log:
        process = subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=environment,
                                   creationflags=subprocess.CREATE_NO_WINDOW)
        try:
            code = process.wait(timeout=min(300,deadline.remaining()-60))
        except subprocess.TimeoutExpired:
            subprocess.run(['taskkill','/PID',str(process.pid),'/T','/F'],capture_output=True)
            raise StopWork('AROSICS stopped at its five-minute execution limit.')
    if code:
        raise StopWork('AROSICS worker failed; inspect its local log, no shift is claimed.')
    result = json.loads((raw/'arosics_result.json').read_text(encoding='utf-8'))
    if not result['accepted_windows']:
        raise StopWork('AROSICS returned no reliable matching windows; no shift is claimed.')
    result['input_mask_rule']=json.loads((raw/'optional_geometry_inputs.json').read_text(encoding='utf-8'))['mask_rule']
    result['input_sha256']={name:file_hash(raw/name) for name in ['coreg_s2_nir.tif','coreg_enmap_nir.tif']}
    write_json(raw/'arosics_result.json',result)
    print(json.dumps(result,indent=2),flush=True)
    return result


def compare_reflectance(cfg, root, deadline):
    import numpy as np
    import pandas as pd
    import rasterio
    from affine import Affine
    from rasterio.warp import reproject, Resampling
    from data_io import read_srf_xlsx, response_weights
    from real_enmap_processing import load_context, warp, refresh_report
    from real_enmap_io import StopWork, write_json, utc_now, file_hash
    raw,meta,prepared,profile,bounds=load_context(cfg,root)
    geometry=json.loads((raw/'arosics_result.json').read_text(encoding='utf-8'))
    audit=json.loads((raw/'s2_original_radiometry_audit.json').read_text(encoding='utf-8'))
    if not geometry['accepted_windows'] or not audit['no_enmap_reference_or_prediction_errors_used']:
        raise StopWork('Geometry or independent S2 radiometry check is missing; no cross-mission accuracy is claimed.')
    if not geometry.get('input_mask_rule','').startswith('Independent'):
        raise StopWork('AROSICS must use independent source masks, not a common artificial nodata pattern.')
    for name,checksum in geometry['input_sha256'].items():
        if file_hash(raw/name)!=checksum:
            raise StopWork('AROSICS input changed after estimation; rerun optional-coreg.')
    calibration=json.loads((raw/'s2_actual_summary.json').read_text(encoding='utf-8'))
    if any(row['effective_offset']!=audit['effective_cog_physical_offset'] for row in calibration['calibration']):
        raise StopWork('S2 preparation changed after its original-product audit; rerun optional-audit.')
    with np.load(raw/'quality_prepared.npz') as z:
        native=z['source_diagnostic_reflectance'];axis=z['source_axis_nm'];native_transform=Affine(*z['native_transform'])
    with np.load(raw/'s2_actual_prepared.npz') as z:
        actual=z['reflectance'];s2_clear=z['clear']
    response_axis,responses=read_srf_xlsx(root/'data/sentinel2_srf.xlsx','S2C',S2_BANDS)
    weights,coverage=response_weights(axis,response_axis,responses,35)
    if np.min(coverage)<.95:
        raise StopWork('Real EnMAP lacks adequate S2C response coverage.')
    shifted_transform=Affine(native_transform.a,native_transform.b,native_transform.c+geometry['median_x_correction_m'],
                             native_transform.d,native_transform.e,native_transform.f+geometry['median_y_correction_m'])
    nominal=np.full_like(actual,np.nan);corrected=nominal.copy();joint=np.zeros_like(actual,dtype=bool)
    rows=[]
    for band in range(len(BANDS)):
        deadline.check(10)
        needed=weights[:,band]>0
        valid=np.isfinite(native[needed])
        supported_weight=np.einsum('bhw,b->hw',valid.astype('float32'),weights[needed,band],optimize=True)
        convolved=np.einsum('bhw,b->hw',np.nan_to_num(native[needed]),weights[needed,band],optimize=True)
        source_good=supported_weight>=1-1e-6
        convolved[~source_good]=np.nan
        good_masks=[]
        for name,transform,output in [('nominal',native_transform,nominal),('arosics',shifted_transform,corrected)]:
            native_profile=profile.copy();native_profile.update(transform=transform,height=native.shape[1],width=native.shape[2])
            good=warp(source_good,native_profile,profile,Resampling.bilinear)>=1-1e-6
            reproject(convolved,output[band],src_transform=transform,src_crs=profile['crs'],
                      dst_transform=profile['transform'],dst_crs=profile['crs'],src_nodata=np.nan,dst_nodata=np.nan,
                      resampling=Resampling.bilinear,num_threads=2)
            good_masks.append(good&np.isfinite(output[band]))
        common=s2_clear&np.isfinite(actual[band])&good_masks[0]&good_masks[1]
        joint[band]=common
        if common.sum()<100:
            raise StopWork('Too few jointly qualified pixels for '+S2_BANDS[band]+'.')
        for geometry_name,array in [('nominal',nominal),('arosics',corrected)]:
            enmap,s2=array[band,common],actual[band,common]
            delta=enmap-s2
            rows.append({'band':S2_BANDS[band],'geometry':geometry_name,'n_pixels':int(common.sum()),
                         'bias_enmap_minus_actual_s2':float(delta.mean()),'rmsd':float(np.sqrt((delta**2).mean())),
                         'pearson_r':float(np.corrcoef(enmap,s2)[0,1]),'enmap_mean':float(enmap.mean()),
                         'actual_s2_mean':float(s2.mean()),'full_s2c_srf_coverage':float(coverage[band])})
    pd.DataFrame(rows).to_csv(raw/'cross_mission_band_metrics.csv',index=False)
    np.savez_compressed(raw/'cross_mission_pairs.npz',enmap_nominal=nominal,enmap_arosics=corrected,actual_s2=actual,joint=joint)
    result={'time_utc':utc_now(),'enmap_date':'2025-05-11','actual_s2_date':'2025-05-12','day_gap':1,
            'actual_s2_platform':'Sentinel-2C','response_platform':'S2C, official ESA SRF v4.0',
            'response_coverage_minimum':float(min(coverage)),
            'geometry_rule':'AROSICS median map-coordinate correction applied only to this comparison; required RF results remain nominal.',
            'paired_mask_rule':'Each band uses the same valid pixels before and after correction; S2 SCL 4/5/6 and source nonzero; EnMAP atmosphere/testflag and every contributing SRF source band good; all bilinear neighbors good.',
            'radiometry_rule':'Original JP2/COG and ESA product-XML audit; no EnMAP-based radiometric fitting.',
            'limitations':'One-day gap, differing atmosphere, view/sun geometry, native PSFs and spatial resampling; not a sensor calibration validation.',
            'metrics':rows}
    write_json(raw/'cross_mission_summary.json',result)
    comparison_figure(cfg,root)
    refresh_report(cfg,root)
    print(pd.DataFrame(rows).to_string(index=False),flush=True)
    return result


def comparison_figure(cfg,root):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    raw=root/cfg['raw_directory']
    metrics=json.loads((raw/'cross_mission_summary.json').read_text(encoding='utf-8'))
    with np.load(raw/'real_predictions.npz') as z:mean=z['enmap_mean'];std=z['enmap_std']
    with np.load(raw/'cross_mission_pairs.npz') as z:enmap=z['enmap_arosics'];s2=z['actual_s2'];good=z['joint']
    fig=plt.figure(figsize=(15,16),constrained_layout=True)
    outer=fig.add_gridspec(4,1,height_ratios=[1.1,1.1,1,1])
    classes=['impervious','vegetation','soil','water']
    for row,arrays,label,limits in [(0,mean,'Frozen RF fraction',(0,1)),(1,std,'Ensemble SD',(0,.15))]:
        inner=outer[row].subgridspec(1,4)
        for i,name in enumerate(classes):
            ax=fig.add_subplot(inner[0,i]);draw=ax.imshow(arrays[i],cmap='viridis',vmin=limits[0],vmax=limits[1])
            ax.set_title(name+': '+label,fontsize=9);ax.set_xticks([]);ax.set_yticks([])
            fig.colorbar(draw,ax=ax,shrink=.7)
    rng=np.random.default_rng(20261104)
    for row in (2,3):
        inner=outer[row].subgridspec(1,5)
        for column in range(5):
            band=(row-2)*5+column
            ax=fig.add_subplot(inner[0,column])
            x,y=s2[band,good[band]],enmap[band,good[band]]
            selected=rng.choice(len(x),min(5000,len(x)),replace=False)
            ax.scatter(x[selected],y[selected],s=2,alpha=.15,color='#376f92',rasterized=True)
            ax.plot([-.05,.8],[-.05,.8],color='#b04b4b',lw=.8)
            ax.set_xlim(-.05,.8);ax.set_ylim(-.05,.8);ax.set_aspect('equal')
            measured=next(m for m in metrics['metrics'] if m['band']==S2_BANDS[band] and m['geometry']=='arosics')
            ax.set_title(f"{S2_BANDS[band]}: n={len(x):,}\nbias={measured['bias_enmap_minus_actual_s2']:.3f}; RMSD={measured['rmsd']:.3f}",fontsize=9)
            ax.set_xlabel('Actual S2C reflectance',fontsize=8)
            if column==0:ax.set_ylabel('Convolved EnMAP reflectance',fontsize=8)
            ax.tick_params(labelsize=7)
    fig.suptitle('Real EnMAP first look — RF maps use strict complete-feature QA (3.09% support)\n'
                 'Scatter: EnMAP 2025-05-11 vs actual S2C 2025-05-12, S2C SRFs, independent per-band QA and AROSICS\n'
                 'Scatter shows at most 5,000 pairs/band and axes -0.05–0.8; statistics use ALL qualified pairs; RF historical comparison has a 16-year gap',fontsize=11)
    fig.savefig(root/cfg['figure_directory']/'transfer_summary.png',dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coreg-worker-raw', required=True)
    coreg_worker(parser.parse_args().coreg_worker_raw)
