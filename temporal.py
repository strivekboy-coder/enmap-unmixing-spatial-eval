"""Cloud-aware real S2 temporal features and independent external-endmember models."""
from __future__ import annotations

import json
import hashlib
import logging
import time
import warnings
import calendar

import dask
import numpy as np
import pandas as pd
import pystac
import rasterio
import xarray as xr
from odc.geo.geobox import GeoBox
from odc.stac import load
from shapely.geometry import box, shape as geometry_shape
from sklearn.ensemble import RandomForestRegressor
from rasterio.warp import transform_bounds

from data_io import write_raster
from evaluation import summarize, block_ids, rho
from mixing import mixtures
from stage_b import dump, stac_search, stop

LOG = logging.getLogger(__name__)
BANDS = ['blue', 'green', 'red', 'rededge1', 'rededge2', 'rededge3', 'nir', 'nir08', 'swir16', 'swir22']


def point_cache_signature(items,rows,cols,tx,shape):
    """Raw cache identity includes actual locations, grid and source assets."""
    payload={'version':1,'shape':list(shape),'transform':list(tx)[:6],
             'assets':[(item['id'],[item['assets'][band]['href'] for band in ['scl']+BANDS]) for item in items],
             'resampling':'scl-nearest;spectra-average','dtype':'uint16',
             'rows_sha256':hashlib.sha256(np.asarray(rows,dtype='int64').tobytes()).hexdigest(),
             'cols_sha256':hashlib.sha256(np.asarray(cols,dtype='int64').tobytes()).hexdigest()}
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()


def calibrate_reflectance(raw, dates, items, bands, audited_zero_offset_items=()):
    """Apply each asset's documented scale/offset, aligned by acquisition date.

    Default: Earth Search raster:bands scale/offset. A small explicit override
    list is permitted only after a same-date red/NIR/SWIR audit proves the
    provider offset was already applied; the original metadata stay in the audit.
    """
    by_day = {item['properties']['datetime'][:10]:item for item in items}
    scales = np.empty((len(dates),len(bands)),dtype='float32')
    offsets = scales.copy()
    audit=[]
    for t,date in enumerate(dates):
        item=by_day[str(date)]
        for b,band in enumerate(bands):
            meta=item['assets'][band].get('raster:bands',[{}])[0]
            scales[t,b]=meta.get('scale',.0001)
            declared=meta.get('offset',0.)
            if item['id'] in audited_zero_offset_items:
                if item['properties'].get('earthsearch:boa_offset_applied') is not True:
                    raise ValueError(f'Audited radiometry override no longer matches provider flag: {item["id"]}')
                offsets[t,b]=0.
            else:
                offsets[t,b]=declared
            audit.append({'date':str(date),'item':item['id'],'band':band,
                          'scale':float(scales[t,b]),'offset':float(offsets[t,b]),
                          'declared_offset':declared,'audited_override':item['id'] in audited_zero_offset_items,
                          'asset_href':item['assets'][band]['href'],
                          'processing_baseline':item['properties'].get('s2:processing_baseline'),
                          'boa_offset_applied_flag':item['properties'].get('earthsearch:boa_offset_applied')})
    return raw.astype('float32')*scales[:,:,None]+offsets[:,:,None],audit


def choose_daily_items(features, bbox):
    """Avoid counting repeat processing products and overlapping swaths as dates."""
    region = box(*bbox)
    groups = {}
    for item in features:
        day = item['properties']['datetime'][:10]
        coverage = geometry_shape(item['geometry']).intersection(region).area / region.area
        cloud = item['properties'].get('eo:cloud_cover', 100.)
        # A primary item covering the tutorial takes precedence over a clipped
        # adjacent swath; then choose lower scene cloud cover deterministically.
        rank = (round(coverage, 4), -cloud, item['id'])
        if day not in groups or rank > groups[day][0]:
            groups[day] = (rank, item)
    return [groups[day][1] for day in sorted(groups)]


def make_cube(items, geobox, bands, cfg):
    chunk = cfg['extension']['chunks_pixels']
    return load([pystac.Item.from_dict(item) for item in items], bands=bands, geobox=geobox,
                chunks={'time': 1, 'x': chunk, 'y': chunk}, groupby='solar_day',
                resampling={b: 'nearest' if b == 'scl' else 'average' for b in bands},
                # odc-loader 0.6.5 passes a dtype mapping to NumPy's scalar
                # chunk resolver. Uniform uint16 avoids that incompatibility;
                # SCL category values and nearest resampling are unchanged.
                dtype='uint16', fail_on_error=True)


def year_observations(year, cfg, out, profile, valid_image):
    bounds = rasterio.transform.array_bounds(profile['height'], profile['width'], profile['transform'])
    bbox = transform_bounds(profile['crs'], 'EPSG:4326', *bounds)
    ext = cfg['extension']
    first_month,last_month=min(ext['months']),max(ext['months'])
    last_day=calendar.monthrange(year,last_month)[1]
    features = stac_search(ext['endpoint'], ext['collection'], bbox,
                           f'{year}-{first_month:02d}-01T00:00:00Z/{year}-{last_month:02d}-{last_day:02d}T23:59:59Z', out/'cube'/f'stac_{year}.json')
    features=[item for item in features if int(item['properties']['datetime'][5:7]) in ext['months']]
    items = choose_daily_items(features, bbox)
    LOG.info('S2 %d: %d STAC products -> %d distinct daily items', year, len(features), len(items))
    if not items:
        return {'year': year, 'products': 0, 'days': 0, 'median_clear_dates': 0., 'p10_clear_dates': 0., 'supported_share': 0.}, items, None, None
    cache = out/'cube'/f'scl_{year}.npz'
    if cache.exists():
        with np.load(cache) as saved:
            values, dates = saved['scl'], saved['dates']
    else:
        geobox = GeoBox((profile['height'], profile['width']), profile['transform'], str(profile['crs']))
        cube = make_cube(items, geobox, ['scl'], cfg)
        LOG.info('Compute %d SCL xarray/Dask cube %s', year, dict(cube.sizes))
        values = cube.scl.compute().values
        dates = cube.time.values.astype('datetime64[D]').astype(str)
        np.savez_compressed(cache, scl=values, dates=dates)
    clear = np.isin(values, ext['clear_scl_values'])
    counts = clear.sum(axis=0)
    observed = values > 0
    count_values = counts[valid_image]
    support = float(np.mean(count_values >= ext['minimum_clear_observations']))
    stats = {'year': year, 'products': len(features), 'days': len(items),
             'median_clear_dates': float(np.median(count_values)), 'p10_clear_dates': float(np.quantile(count_values, .1)),
             'minimum_clear_dates': int(count_values.min()), 'maximum_clear_dates': int(count_values.max()),
             'supported_share': support}
    clear_fraction = clear[:, valid_image].mean(axis=1)
    observed_fraction = observed[:, valid_image].mean(axis=1)
    daily = pd.DataFrame({'year': year, 'date': dates, 'clear_area_fraction': clear_fraction,
                          'observed_area_fraction': observed_fraction})
    by_day={item['properties']['datetime'][:10]:item for item in items}
    daily['selected_item']=[by_day[str(date)]['id'] for date in dates]
    daily['scene_cloud_cover_percent']=[by_day[str(date)]['properties'].get('eo:cloud_cover',np.nan) for date in dates]
    daily['usable'] = ((daily.clear_area_fraction >= ext['usable_clear_area_fraction']) &
                        (daily.observed_area_fraction >= ext['minimum_observed_area_fraction']))
    daily['month'] = pd.to_datetime(daily.date).dt.month
    daily.to_csv(out/f'daily_cloud_{year}.csv', index=False)
    LOG.info('S2 %d clear support: %s', year, stats)
    return stats, items, counts, daily


def acquisition_outputs(cfg, out, profile, valid_image, selected, observations):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    records = []
    for stats, items, counts, daily in observations:
        if daily is None:
            continue
        for month in cfg['extension']['months']:
            rows = daily[daily.month == month]
            p = float(rows.usable.mean()) if len(rows) else np.nan
            scene_values=rows.scene_cloud_cover_percent.dropna()
            scene_p=float((scene_values<=cfg['extension']['usable_scene_cloud_percent']).mean()) if len(scene_values) else np.nan
            records.append({'year': stats['year'], 'month': month, 'distinct_dates': len(rows),
                            'usable_dates': int(rows.usable.sum()), 'local_usable_probability': p,
                            'expected_enmap_attempts_proxy': 1/p if p > 0 else np.inf,
                            'scene_low_cloud_probability':scene_p,
                            'scene_expected_attempts_proxy':1/scene_p if scene_p>0 else np.inf,
                            'mean_local_clear_area_fraction': rows.clear_area_fraction.mean() if len(rows) else np.nan})
    monthly = pd.DataFrame(records)
    monthly.to_csv(out/'monthly_acquisition_proxy.csv', index=False)
    chosen = next(result for result in observations if result[0]['year'] == selected)
    stats, _, counts, daily = chosen
    seasonal_p = float(daily.usable.mean())
    dump(out/'seasonal_acquisition_proxy.json', {'year': selected, 'usable_probability': seasonal_p,
         'expected_attempts': 1/seasonal_p if seasonal_p else None,
         'assumptions': 'S2 local SCL cloudiness proxies EnMAP; acquisition attempts are independent with stationary success probability. Revisit, geometry, seasonal persistence and scheduling constraints are ignored.'})
    write_raster(out/f'clear_observation_count_{selected}.tif', np.where(valid_image, counts, np.nan)[None], profile, ['clear_distinct_dates'])
    extent = [0, profile['width']*.03, 0, profile['height']*.03]
    fig, axes = plt.subplots(1, 2, figsize=(11, 6), layout='constrained')
    im = axes[0].imshow(np.where(valid_image, counts, np.nan), extent=extent, cmap='viridis')
    axes[0].set_title(f'{selected}: clear SCL observations\nmedian {stats["median_clear_dates"]:.0f} dates / pixel')
    axes[0].set_xlabel('Easting from grid origin (km)'); axes[0].set_ylabel('Northing from southern edge (km)')
    fig.colorbar(im, ax=axes[0], label='Distinct clear dates', shrink=.7)
    months = cfg['extension']['months']
    years = cfg['extension']['candidate_years']
    for i, year in enumerate(years):
        rows = monthly[monthly.year == year].set_index('month').reindex(months)
        axes[1].bar(np.arange(len(months)) + (i-.5)*.35, rows.local_usable_probability, width=.35, label=str(year))
    axes[1].set_xticks(np.arange(len(months)), ['Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct'])
    axes[1].set_ylim(0, 1); axes[1].set_ylabel('Probability of >=80% local clear area')
    axes[1].set_title('Local clear-acquisition proxy\nAt least 95% observed footprint required')
    axes[1].legend()
    fig.suptitle('Acquisition planning from real S2 cloud masks | proxy, not an EnMAP scheduling forecast')
    fig.savefig(out/'clear_observations_acquisition.png', dpi=170, bbox_inches='tight'); plt.close(fig)
    return monthly


def load_point_features(cfg, root, out, selected, items, candidates, tx, shape, valid_image, profile):
    """Dask loads only chunks intersecting evaluation pixels or library points."""
    pad = int(round((profile['transform'].c-tx.c)/30))
    er, ec = np.where(valid_image)
    rows = np.r_[er+pad, candidates.row.to_numpy()]
    cols = np.r_[ec+pad, candidates.col.to_numpy()]
    points = xr.DataArray(np.arange(len(rows)), dims='point')
    point_rows = xr.DataArray(rows, dims='point')
    point_cols = xr.DataArray(cols, dims='point')
    n_eval = len(er)
    path = out/'cube'/f'spectral_points_{selected}.npz'
    provenance_path=out/'cube'/f'spectral_points_{selected}_provenance.json'
    signature=point_cache_signature(items,rows,cols,tx,shape)
    cached_signature=json.loads(provenance_path.read_text(encoding='utf-8')).get('raw_cache_signature') if provenance_path.exists() else None
    if path.exists() and cached_signature==signature:
        with np.load(path) as data:
            spectral, scl, dates = data['spectral'], data['scl'], data['dates']
    else:
        cube = make_cube(items, GeoBox(shape, tx, str(profile['crs'])), ['scl']+BANDS, cfg)
        dates = cube.time.values.astype('datetime64[D]').astype(str)
        LOG.info('Load SCL at %d image and %d external points', n_eval, len(candidates))
        scl = cube.scl.isel(y=point_rows, x=point_cols).transpose('time', 'point').compute().values
        spectral = np.empty((len(dates), len(BANDS), len(rows)), dtype='uint16')
        for index, band in enumerate(BANDS):
            tic = time.perf_counter()
            LOG.info('Load real S2 band %s (%d/%d) at supported points', band, index+1, len(BANDS))
            spectral[:, index] = cube[band].isel(y=point_rows, x=point_cols).transpose('time', 'point').compute().values
            LOG.info('Loaded %s in %.1fs', band, time.perf_counter()-tic)
        partial=path.with_name(path.stem+'.partial.npz')
        np.savez_compressed(partial, spectral=spectral, scl=scl, dates=dates)
        partial.replace(path)
        dump(provenance_path,{'raw_cache_signature':signature,'image_points':n_eval,'external_points':len(candidates),
                              'note':'Locations, source asset URLs, grid and resampling are part of the cache key.'})
    overrides=cfg['extension'].get('audited_offset_already_applied_items',[])
    if overrides:
        audit=pd.read_csv(out/'same_date_radiometric_audit.csv')
        for item_id in overrides:
            rows=audit[audit.selected_item==item_id]
            if (not {'red','nir','swir16'}.issubset(set(rows.band)) or rows.n_pixels.min()<1000 or
                 rows.median_raw_DN_difference.abs().max()>cfg['extension']['maximum_same_date_median_DN_difference']):
                raise ValueError(f'Radiometric offset override lacks a passing independent audit: {item_id}')
    scaled,calibration=calibrate_reflectance(spectral,dates,items,BANDS,overrides)
    pd.DataFrame(calibration).to_csv(out/'radiometric_calibration_by_date_band.csv',index=False)
    unique=sorted({(row['scale'],row['offset']) for row in calibration})
    dump(out/'radiometric_calibration.json', {'unique_scale_offset':unique,'bands':BANDS,
         'source':'https://github.com/Element84/earth-search/blob/main/docs/collections/sentinel-2-l2a.md#gainoffset-in-items-after-jan-25-2022',
         'audited_offset_already_applied_items':overrides,
         'note':'Apply asset scale/offset per date, except five explicitly audited items whose -0.1 metadata would double-apply the already subtracted offset. Same-date red/NIR/SWIR comparisons use no reference labels.'})
    clear = np.isin(scl, cfg['extension']['clear_scl_values'])
    good = clear & (spectral > 0).all(axis=1)
    external_coverage = good[:, n_eval:].mean(axis=1)
    best_index = int(external_coverage.argmax())
    scaled = np.where(good[:, None, :], scaled, np.nan)
    red, nir = scaled[:, 2], scaled[:, 6]
    ndvi = np.where((red>=0)&(nir>=0)&(nir+red>1e-6),(nir-red)/(nir+red+1e-8),np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        medians = np.nanmedian(scaled, axis=0).T
        ndvi_features = np.array([np.nanmedian(ndvi, axis=0), np.nanmax(ndvi, axis=0),
                                  np.nanmax(ndvi, axis=0)-np.nanmin(ndvi, axis=0)]).T
    temporal = np.column_stack([medians, ndvi_features])
    # Identical feature schema for one-date and seasonal views: with one date,
    # each band median is that observation, NDVI median=max, and amplitude=0.
    # This avoids giving only the temporal model an engineered NDVI ratio.
    one_ndvi=ndvi[best_index]
    single=np.column_stack([scaled[best_index].T,one_ndvi,one_ndvi,
                            np.where(np.isfinite(one_ndvi),0.,np.nan)])
    sufficient = good.sum(axis=0) >= cfg['extension']['minimum_clear_observations']
    temporal[~sufficient]=np.nan
    common = sufficient & np.isfinite(temporal).all(axis=1) & np.isfinite(single).all(axis=1)
    metadata = {'selected_year': selected, 'best_date': dates[best_index], 'best_date_external_clear_share': float(external_coverage[best_index]),
                'feature_names': [band+'_seasonal_median' for band in BANDS]+['ndvi_median', 'ndvi_max', 'ndvi_amplitude'],
                'single_feature_names':[band+'_best_date' for band in BANDS]+['ndvi_one_date_median','ndvi_one_date_max','ndvi_one_date_amplitude'],
                'single_date_schema':'Same 13 features as temporal model; one-date NDVI median=max and amplitude=0.',
                'n_eval_image_points': n_eval, 'n_external_candidates': len(candidates),
                'n_external_common_supported': int(common[n_eval:].sum()),
                'n_image_common_supported': int(common[:n_eval].sum())}
    dump(out/'temporal_feature_manifest.json', metadata)
    return single, temporal, common, good.sum(axis=0), n_eval, metadata


def temporal_analysis(cfg, root, out, profile, reference, valid, evidence, stable, candidates, tx, shape):
    with np.load(root/cfg['results']/'cache/prepared.npz') as data:
        valid_image = data['valid_image']
    with dask.config.set(scheduler='threads', num_workers=cfg['extension']['dask_workers']):
        observations = [year_observations(year, cfg, out, profile, valid_image) for year in cfg['extension']['candidate_years']]
        pd.DataFrame([r[0] for r in observations]).to_csv(out/'year_selection.csv', index=False)
        chosen = max(observations, key=lambda r: (r[0]['median_clear_dates'], r[0]['p10_clear_dates'], -r[0]['year']))
        selected = chosen[0]['year']
        acquisition_outputs(cfg, out, profile, valid_image, selected, observations)
        evidence['selected_s2_year'] = selected
        evidence['s2_minus_2009_years'] = selected-2009
        evidence['worldcover_minus_s2_years'] = cfg['extension']['worldcover_year']-selected
        evidence['clear_observation_stats'] = chosen[0]
        if chosen[0]['supported_share'] < cfg['extension']['minimum_supported_pixel_share']:
            return stop(out, f"Only {chosen[0]['supported_share']:.1%} of image pixels have >=5 clear dates, below the fixed 80% gate.", evidence)
        single, temporal, common, counts, n_eval, feature_meta = load_point_features(
            cfg, root, out, selected, chosen[1], candidates, tx, shape, valid_image, profile)
    library_ok = common[n_eval:]
    supported = candidates[library_ok].copy()
    supported['clear_observation_count'] = counts[n_eval:][library_ok]
    supported.to_csv(out/'external_endmembers_used.csv', index=False)
    support = supported.groupby('class').size().reindex(cfg['classes'], fill_value=0)
    evidence['endmembers_after_feature_filter'] = support.to_dict()
    dump(out/'preflight_all_gates.json', evidence)
    if support.min() < cfg['extension']['minimum_endmembers_per_class']:
        return stop(out, f'External endmembers after clear-date/best-date filtering fail the >=30/class gate: {support.to_dict()}', evidence)
    from temporal_models import train_temporal_models
    return train_temporal_models(cfg, root, out, profile, reference, valid, valid_image, stable, supported,
                                  single, temporal, common, n_eval, evidence, feature_meta)
