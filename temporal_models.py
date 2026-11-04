"""Same-source temporal comparison and fixed, untrained late fusion."""
from __future__ import annotations

import logging
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sklearn.ensemble import RandomForestRegressor

from data_io import write_raster
from evaluation import summarize
from mixing import mixtures
from plotting import block_grid
from report import markdown_table
from stage_b import dump, stop

LOG = logging.getLogger(__name__)


def train_temporal_models(cfg, root, out, profile, reference, valid, valid_image, stable, supported,
                          single, temporal, common, n_eval, evidence, feature_meta):
    n = cfg['mixing']['ensemble_size']
    library_indices = supported.index.to_numpy()+n_eval
    labels = supported.class_index.to_numpy()
    joined = np.column_stack([single[library_indices], temporal[library_indices]])
    feature_names = feature_meta['single_feature_names']+feature_meta['feature_names']
    feature_table = supported.reset_index(drop=True).copy()
    for col, name in enumerate(feature_names):
        feature_table[name] = joined[:, col]
    feature_table.to_csv(out/'external_endmember_features.csv', index=False)
    pred = {name: np.full((n, 4, n_eval), np.nan, dtype='float32') for name in ['c0_S2_best_date', 'c_S2_temporal']}
    features = {'c0_S2_best_date': single[:n_eval], 'c_S2_temporal': temporal[:n_eval]}
    valid_features = {name: np.isfinite(values).all(axis=1) for name, values in features.items()}
    timings = []
    for run in range(n):
        for target, name in enumerate(cfg['classes']):
            seed = cfg['seed']+run*100+target
            x, y = mixtures(joined, labels, target, cfg['mixing'], seed)
            for model_name in pred:
                tic = time.perf_counter()
                train = x[:, :single.shape[1]] if model_name == 'c0_S2_best_date' else x[:, single.shape[1]:]
                model = RandomForestRegressor(**cfg['forest'], random_state=seed)
                model.fit(train, y)
                mask = valid_features[model_name]
                pred[model_name][run, target, mask] = model.predict(features[model_name][mask])
                timings.append({'model': model_name, 'run': run+1, 'class': name,
                                'seconds': time.perf_counter()-tic, 'synthetic_oob_r2': model.oob_score_, 'seed': seed})
                LOG.info('Temporal %s run %d class %s: %.2fs', model_name, run+1, name, timings[-1]['seconds'])
    with np.load(root/cfg['results']/'cache/ensemble_predictions.npz') as stage_a:
        pred['a_EnMAP_library'] = stage_a['EnMAP'][:n]
        pred['b_S2_like_library'] = stage_a['S2-like'][:n]
    pred['d_late_fusion'] = (pred['a_EnMAP_library']+pred['c_S2_temporal'])/2
    np.savez_compressed(out/'cube/temporal_model_members.npz', **pred)
    pd.DataFrame(timings).to_csv(out/'model_runtime.csv', index=False)
    means, stds = {}, {}
    for name, values in pred.items():
        mean = np.full((4, profile['height'], profile['width']), np.nan, dtype='float32')
        std = mean.copy()
        mean[:, valid_image] = values.mean(axis=0)
        std[:, valid_image] = values.std(axis=0, ddof=0)
        means[name], stds[name] = mean, std
        write_raster(out/f'{name}_fractions.tif', mean, profile, cfg['classes'])
        write_raster(out/f'{name}_uncertainty.tif', std, profile, cfg['classes'])
    supported_image = np.zeros(valid.shape, dtype=bool)
    supported_image[valid_image] = common[:n_eval]
    full = valid & supported_image
    screened = full & stable
    if not full.any() or not screened.any():
        return stop(out, 'No common supported reference pixels remain after fixed feature/consistency masks.', evidence)
    write_raster(out/'evaluation_common_support_mask.tif',np.where(valid,full,np.nan)[None],profile,['common_supported_reference'])
    write_raster(out/'evaluation_consistency_support_mask.tif',np.where(valid,screened,np.nan)[None],profile,['consistent_supported_reference'])
    reference_dominant=reference.argmax(axis=0)
    support_rows=[]
    for target,name in enumerate(cfg['classes']):
        original=valid&(reference_dominant==target)
        support_rows.append({'dominant_reference':name,'original_pixels':int(original.sum()),
                              'common_supported_pixels':int((original&full).sum()),
                              'consistent_supported_pixels':int((original&screened).sum())})
    pd.DataFrame(support_rows).to_csv(out/'evaluation_support_by_context.csv',index=False)
    all_summary, all_blocks, all_contexts = [], [], []
    for subset, mask in [('full_reference_supported', full), ('worldcover_consistent_supported', screened)]:
        summary, blocks, contexts = summarize(cfg, reference, mask, means, stds, n)
        for table in [summary, blocks, contexts]:
            table.insert(0, 'evaluation_subset', subset)
        all_summary.append(summary); all_blocks.append(blocks); all_contexts.append(contexts)
    summary = pd.concat(all_summary, ignore_index=True)
    blocks = pd.concat(all_blocks, ignore_index=True)
    contexts = pd.concat(all_contexts, ignore_index=True)
    summary.to_csv(out/'summary.csv', index=False)
    blocks.to_csv(out/'blocks.csv', index=False)
    contexts.to_csv(out/'contexts.csv', index=False)
    evidence['reference_support'] = {'original_reference_pixels': int(valid.sum()), 'common_supported_pixels': int(full.sum()),
                                      'consistency_screen_supported_pixels': int(screened.sum()),
                                      'cloud_feature_exclusions': int(valid.sum()-full.sum())}
    pairs = [('c0_S2_best_date','c_S2_temporal'), ('c_S2_temporal','d_late_fusion'), ('a_EnMAP_library','d_late_fusion')]
    tables = []
    for baseline, candidate in pairs:
        a = blocks[blocks.sensor == baseline]
        b = blocks[blocks.sensor == candidate]
        table = a.merge(b, on=['evaluation_subset','class','block_id'], suffixes=('_baseline','_candidate'))
        table['baseline'] = baseline; table['candidate'] = candidate
        table['mae_reduction'] = table.mae_baseline-table.mae_candidate
        tables.append(table)
    gains = pd.concat(tables, ignore_index=True)
    gains.to_csv(out/'value_of_acquisition_blocks.csv', index=False)
    stage_b_figures(cfg, out, profile, means, gains)
    metric = summary[summary.evaluation_subset == 'full_reference_supported'].pivot(index='class', columns='sensor', values='mae').reindex(cfg['classes']+['macro'])
    consistent_metric = summary[summary.evaluation_subset == 'worldcover_consistent_supported'].pivot(index='class', columns='sensor', values='mae').reindex(cfg['classes']+['macro'])
    metric.reset_index().to_csv(out/'mae_model_comparison.csv', index=False)
    year = evidence['selected_s2_year']
    library_counts=pd.DataFrame([{'class':name,'worldcover_year':cfg['extension']['worldcover_year'],
                                  'used_endmembers':int(evidence['endmembers_after_feature_filter'][name])}
                                 for name in cfg['classes']])
    note = f'''# Stage B — method demonstration under a known time gap

**Dates: HyMap/simulated EnMAP 2009; real S2 growing season {year} (gap {year-2009} years); external WorldCover labels {cfg['extension']['worldcover_year']} (another {cfg['extension']['worldcover_year']-year} years after S2). These results are not a same-date sensor comparison.**

The fixed support/cloud/agreement gates passed. The external library lies outside the entire tutorial extent plus a 1 km exclusion, within a 10 km buffer. WorldCover dominant-class agreement is {evidence['worldcover_reference_agreement']:.1%}: {evidence['consistent_pixels']}/{evidence['reference_pixels']} pixels before cloud-driven exclusions. This consistency screen is not proof of temporal stability and preferentially selects agreement with the map used to label the training library; the full-area result below is therefore primary.

External endmembers retained after fixed cloud/feature filtering:

{markdown_table(library_counts)}

## Primary matched-support results

All five models use the same **{int(full.sum())}** qualifying reference pixels, with {int(valid.sum()-full.sum())} excluded because of missing/cloud-limited single-date/seasonal support. Fractions and MAE are on a 0–1 scale. Fixed Stage A models are re-evaluated on this common support, not substituted with the original full-support Stage A metric.

`evaluation_support_by_context.csv` gives the retained counts by dominant class, and the two evaluation-mask GeoTIFFs make all exclusions inspectable.

{markdown_table(metric.reset_index())}

## WorldCover-consistency sensitivity

{int(screened.sum())} common-supported pixels remain under the predeclared WorldCover/reference class-agreement mask.

{markdown_table(consistent_metric.reset_index())}

## Method and provenance

- Years 2017 and 2018 were both queried over Apr–Oct, deduplicated to one completely covering/lower-cloud item per date and compared by median then 10th-percentile clear SCL counts; selection did not use reference errors. `year_selection.csv`, `daily_cloud_*.csv` and `cube/stac_*.json` record this.
- STAC + odc-stac produced an xarray/Dask cube on the exact 30 m EnMAP grid. SCL 4/5/6 were retained with nearest-neighbour resampling; spectral bands used average resampling. Metadata scale/offset are checked before feature extraction. `radiometric_calibration.json` records them.
- External candidate labels use [ESA WorldCover](https://planetarycomputer.microsoft.com/dataset/esa-worldcover), CC BY 4.0. The 10 m labels were sampled on an aligned 10 m UTM grid and aggregated 3×3. At >=95% mapped purity, all nine samples must agree. Class crosswalk: built-up -> impervious; bare/sparse -> soil; permanent water -> water; trees/shrubs/grass/crops/wetlands/moss -> vegetation. Built-up is not identical to impervious cover, and mapped purity is not physical spectral purity.
- (c0) uses the single date **{feature_meta['best_date']}**, selected by external-library clear coverage only. Both (c0)/(c) use the same 13-feature schema: ten band medians plus NDVI median, maximum and amplitude. For the one-date view, each band median is that date's value, NDVI median=max and amplitude=0; the temporal view summarizes the season. Thus NDVI feature engineering is available to both models. They use exactly the same external endmember locations after common cloud/feature filtering. `external_endmembers_used.csv` and `external_endmember_features.csv` list counts, coordinates and actual features. Neither model trains on tutorial reference fractions.
- Synthetic temporal mixing operates in this summary-feature space as specified in the design. NDVI, medians and amplitude are nonlinear functions of reflectance; linearly mixing those features is a modelling approximation, not exact physical mixing of complete multitemporal reflectance sequences.
- Target-first mixing, original endmembers, complexity probabilities and RF settings match Stage A. Three synthetic-training/model seeds define ensemble spread. (d) averages corresponding members of (a) and (c) with fixed 0.5/0.5 weights and no additional training. This is a direct late-fusion sensitivity, not an optimized combination.
- Original library (a)/(b) and external-library (c0)/(c) have different source spectra and class counts. The clean temporal comparison is **c0 vs c**, and the fixed addition-of-hyperspectral comparison is **c vs d**; other contrasts also include library/domain changes. Historical class changes and WorldCover label errors can make an apparently pure external label wrong for the selected S2 year.
- `summary.csv`, `blocks.csv` and `contexts.csv` report per-class metrics, dominant-reference/impervious-gradient contexts and pixel/block spread–error rank correlations for both masks. No independent-pixel significance is claimed. The `value_of_acquisition_blocks.csv` differences and map require reference truth and are retrospective.
- Adding the **2009** simulated hyperspectral map to **{year}** S2 is confounded by their different dates and the historical reference. A reduction in old-reference error is not a measured benefit of acquiring a new hyperspectral scene in the S2 year or today.
- Five reprocessed assets advertised a -0.1 offset despite already offset-corrected pixel values. Same-date legacy comparisons over 71,158 pixels in red/NIR/SWIR confirmed this independently of reference labels. Only the explicitly audited item IDs use effective offset zero; both declared and effective offsets appear in `radiometric_calibration_by_date_band.csv`. `same_date_radiometric_audit.csv` and the fixed QC amendment document the evidence.

## Acquisition planning

`monthly_acquisition_proxy.csv` reports distinct dates, the probability that >=80% of the local valid footprint is clear and >=95% is observed, and `1/p` expected attempts. S2 cloudiness proxies EnMAP: independent attempts, stationary success probability, matched geographic conditions and unconstrained scheduling are assumptions, not established mission performance. A zero-success month has no finite plug-in attempt estimate. Same-day products are not treated as independent acquisition opportunities.

Clear dates here mean strict SCL classes 4/5/6. The spatial counts also reflect SCL classification/shadow behaviour and may depend on surface type; they are not a pure meteorological cloud record. A second, scene-level <=20% cloud-cover proxy is reported alongside the local SCL proxy, rather than treating the two as interchangeable.

![Clear dates and acquisition proxy](clear_observations_acquisition.png)
![Temporal and fusion gains](value_of_acquisition.png)
![Temporal and fused fractions](temporal_fractions.png)

Stage C scaling, official QGIS comparison and real EnMAP transfer were not run.
'''
    (out/'STATUS.md').write_text(note, encoding='utf-8')
    dump(out/'gate_status.json', {'status': 'complete', 'evidence': evidence,
                                 'framing': 'Cross-year method demonstration, not a same-date accuracy comparison.'})
    from finalize_stage_b import finalize
    finalize(cfg,root,root/cfg['results'])
    LOG.info('Stage B complete: %d common reference pixels; %d screened pixels', full.sum(), screened.sum())
    return {'status': 'complete', 'evidence': evidence}


def stage_b_figures(cfg, out, profile, means, gains):
    h, w = profile['height'], profile['width']
    extent = [0, w*.03, 0, h*.03]
    rows = gains[(gains.evaluation_subset == 'full_reference_supported') & (gains['class'] == 'macro')]
    limit = max(.01, rows.mae_reduction.abs().max()*100)
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    fig, axes = plt.subplots(1, 3, figsize=(11, 6), layout='constrained')
    labels = ['Temporal value: MAE(c0) - MAE(c)', 'Added HSI: MAE(c) - MAE(d)', 'Added time series: MAE(a) - MAE(d)']
    for ax, ((baseline, candidate), title) in zip(axes, zip(
        [('c0_S2_best_date','c_S2_temporal'), ('c_S2_temporal','d_late_fusion'), ('a_EnMAP_library','d_late_fusion')], labels)):
        table = rows[(rows.baseline == baseline) & (rows.candidate == candidate)]
        grid = block_grid(table, 'mae_reduction', (h,w), cfg['evaluation']['block_size_pixels'])
        im = ax.imshow(grid, extent=extent, cmap='RdBu', norm=norm, interpolation='nearest')
        ax.set_title(title, fontsize=10); ax.set_facecolor('#eeeeee'); ax.set_xlabel('Easting (km)')
    axes[0].set_ylabel('Northing (km)')
    fig.colorbar(im, ax=axes.tolist(), label='Macro-MAE reduction (percentage points)', shrink=.7)
    fig.suptitle('Retrospective acquisition value | blue helps; red hurts\nKnown 2009 / S2-year / 2021 label gap; reference truth required')
    fig.savefig(out/'value_of_acquisition.png', dpi=170, bbox_inches='tight'); plt.close(fig)
    fig, axes = plt.subplots(2,4,figsize=(12,11),layout='constrained',sharex=True,sharey=True)
    for row, model in enumerate(['c_S2_temporal','d_late_fusion']):
        for target, name in enumerate(cfg['classes']):
            im = axes[row,target].imshow(means[model][target],extent=extent,cmap='viridis',vmin=0,vmax=1)
            axes[row,target].set_title(name if row == 0 else '')
            if target == 0:
                axes[row,target].set_ylabel(model+'\nNorthing (km)')
            if row == 1:
                axes[row,target].set_xlabel('Easting (km)')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Fraction',shrink=.6)
    fig.suptitle('External-library real S2 temporal fractions and fixed late fusion\nCross-year method demonstration')
    fig.savefig(out/'temporal_fractions.png',dpi=170,bbox_inches='tight');plt.close(fig)
