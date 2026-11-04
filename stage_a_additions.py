"""Frozen-baseline spatial gains, water diagnosis, and synthetic-only SVR sensitivity."""
from __future__ import annotations

import json
import logging
import sqlite3
import time
from itertools import product

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import rasterio
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

from data_io import write_raster, sha256
from evaluation import summarize, differences, block_ids, rho
from mixing import mixtures
from report import markdown_table

LOG = logging.getLogger(__name__)


def maps_from_members(members, valid_image, shape, size):
    means, stds = {}, {}
    for sensor, values in members.items():
        mean = np.full((4, *shape), np.nan, dtype='float32')
        std = mean.copy()
        mean[:, valid_image] = values[:size].mean(axis=0)
        std[:, valid_image] = values[:size].std(axis=0, ddof=0)
        means[sensor], stds[sensor] = mean, std
    return means, stds


def gain_diagnostics(cfg, out, results, reference, valid, means):
    size = cfg['mixing']['ensemble_size']
    ids = block_ids(valid.shape, cfg['evaluation']['block_size_pixels'])[valid]
    a = np.abs(means['EnMAP'][:, valid]-reference[:, valid])
    b = np.abs(means['S2-like'][:, valid]-reference[:, valid])
    block_rows = []
    for target, name in list(enumerate(cfg['classes']))+[(None,'macro')]:
        delta = (b-a).mean(axis=0) if target is None else b[target]-a[target]
        for bid in np.unique(ids):
            mask = ids == bid
            block_rows.append({'class': name, 'block_id': int(bid), 'n_pixels': int(mask.sum()),
                               's2_minus_enmap_mae': float(delta[mask].mean()),
                               'enmap_better': bool(delta[mask].mean() > 0),
                               'passes_block_support': bool(mask.sum() >= cfg['evaluation']['minimum_pixels_per_block'])})
    blocks = pd.DataFrame(block_rows)
    blocks.to_csv(out/'rf_block_gains_all_support.csv', index=False)
    ci = pd.read_csv(results/'paired_block_bootstrap.csv')
    ci = ci[ci.ensemble_size == size].copy()
    for index, row in ci.iterrows():
        selected = blocks[blocks['class'] == row['class']]
        supported = selected[selected.passes_block_support]
        ci.loc[index, 'share_all_blocks_enmap_better'] = selected.enmap_better.mean()
        ci.loc[index, 'share_supported_blocks_enmap_better'] = supported.enmap_better.mean()
        ci.loc[index, 'supported_blocks'] = len(supported)
        ci.loc[index, 'bootstrap_resamples'] = cfg['evaluation']['bootstrap_resamples']
    ci.to_csv(out/'rf_gain_intervals_and_block_shares.csv', index=False)
    contexts = pd.read_csv(results/'contexts.csv')
    contexts = contexts[contexts.ensemble_size == size]
    keys = ['class', 'context_type', 'context']
    ca = contexts[contexts.sensor == 'EnMAP'][keys+['n_pixels','mae']]
    cb = contexts[contexts.sensor == 'S2-like'][keys+['n_pixels','mae']]
    context_gain = ca.merge(cb, on=keys, suffixes=('_enmap','_s2'), validate='one_to_one')
    context_gain['s2_minus_enmap_mae'] = context_gain.mae_s2-context_gain.mae_enmap
    extra = []
    edges = [0,.25,.5,.75,1.000001]
    vegetation = reference[1,valid]
    for group, (lo, hi) in enumerate(zip(edges[:-1],edges[1:])):
        mask = (vegetation >= lo) & (vegetation < hi)
        if not mask.any():
            continue
        for target, name in list(enumerate(cfg['classes']))+[(None,'macro')]:
            aa = a.mean(axis=0) if target is None else a[target]
            bb = b.mean(axis=0) if target is None else b[target]
            extra.append({'class':name,'context_type':'reference_vegetation_fraction',
                           'context':['0-25%','25-50%','50-75%','75-100%'][group],
                           'n_pixels_enmap':int(mask.sum()),'mae_enmap':aa[mask].mean(),
                           'n_pixels_s2':int(mask.sum()),'mae_s2':bb[mask].mean(),
                           's2_minus_enmap_mae':(bb[mask]-aa[mask]).mean()})
    context_gain = pd.concat([context_gain,pd.DataFrame(extra)],ignore_index=True)
    context_gain.to_csv(out/'rf_gain_by_context.csv',index=False)
    selected = context_gain[context_gain.context_type == 'dominant_reference']
    matrix = selected.pivot(index='context',columns='class',values='s2_minus_enmap_mae').reindex(index=cfg['classes'],columns=cfg['classes']+['macro'])*100
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    limit=max(.01,np.abs(matrix.to_numpy()).max())
    im=ax.imshow(matrix,cmap='RdBu',norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit),aspect='auto')
    ax.set_xticks(range(5),matrix.columns);ax.set_yticks(range(4),matrix.index)
    ax.set_xlabel('Predicted fraction');ax.set_ylabel('Dominant reference class')
    for r in range(4):
        for c in range(5):
            value=matrix.iloc[r,c]
            ax.text(c,r,f'{value:+.2f}',ha='center',va='center',color='white' if abs(value)>.65*limit else 'black')
    fig.colorbar(im,ax=ax,label='S2-like MAE - EnMAP MAE (pp)')
    ax.set_title('Gain by land-cover context | blue favours EnMAP\nSparse soil/water contexts are diagnostic, not independent replications')
    fig.savefig(out/'gain_by_context.png',dpi=180,bbox_inches='tight');plt.close(fig)
    return ci,context_gain


def water_diagnosis(cfg,out,root,image,library,labels,wavelengths,valid_image,reference,valid,means,profile):
    settings=cfg['stage_a_additions']
    scale=settings['diagnostic_input_scale']
    vnir=np.full(valid.shape,np.nan,dtype='float32')
    swir=vnir.copy()
    vnir[valid_image]=image[wavelengths<1000][:,valid_image].mean(axis=0)/scale
    swir[valid_image]=image[wavelengths>1500][:,valid_image].mean(axis=0)/scale
    q=settings['dark_surface_quantile']
    thresholds=[float(np.quantile(values[valid_image],q)) for values in [vnir,swir]]
    shadow_like=(vnir<thresholds[0])&(swir<thresholds[1])&valid_image
    brightness=(vnir+swir)/2
    quartiles=np.quantile(brightness[valid_image],[.25,.5,.75])
    brightness_group=np.digitize(brightness[valid],quartiles)
    refwater=reference[3,valid]
    dominant=reference[:,valid].argmax(axis=0)
    rows,correlations=[],[]
    for sensor,pred in means.items():
        water=pred[3,valid]
        signed=water-refwater
        false=(refwater<=settings['water_absence_threshold'])&(water>=settings['water_false_prediction_threshold'])
        for name,b in [('mean_VNIR',vnir[valid]),('mean_SWIR',swir[valid])]:
            for subset,mask in [('all',np.ones(len(water),bool)),('reference_water_absent',refwater<=settings['water_absence_threshold'])]:
                correlations.append({'sensor':sensor,'brightness':name,'subset':subset,'n_pixels':int(mask.sum()),
                                      'spearman_signed_water_error':rho(b[mask],signed[mask]),
                                      'spearman_absolute_water_error':rho(b[mask],np.abs(signed[mask])),
                                      'spearman_predicted_water':rho(b[mask],water[mask])})
        groups=[('dominant_reference',dominant,cfg['classes']),
                 ('shadow_like_proxy',shadow_like[valid].astype(int),['other','dark_VNIR_and_SWIR']),
                 ('image_brightness_quartile',brightness_group,['Q1_darkest','Q2','Q3','Q4_brightest'])]
        for kind,values,names in groups:
            for code,name in enumerate(names):
                mask=values==code
                if not mask.any():continue
                rows.append({'sensor':sensor,'context_type':kind,'context':name,'n_pixels':int(mask.sum()),
                             'mean_predicted_water':float(water[mask].mean()),'mean_reference_water':float(refwater[mask].mean()),
                             'water_bias':float(signed[mask].mean()),'water_mae':float(np.abs(signed[mask]).mean()),
                             'false_water_pixels':int(false[mask].sum()),'false_water_pixel_share':float(false[mask].mean()),
                             'mean_VNIR':float(vnir[valid][mask].mean()),'mean_SWIR':float(swir[valid][mask].mean())})
    groups=pd.DataFrame(rows);corr=pd.DataFrame(correlations)
    groups.to_csv(out/'water_diagnosis_contexts.csv',index=False)
    corr.to_csv(out/'water_brightness_correlations.csv',index=False)
    _,water_targets=mixtures(library,labels,3,cfg['mixing'],cfg['seed']+3)
    (out/'water_synthetic_target_distribution.json').write_text(json.dumps({
        'n_training_rows':len(water_targets),'zero_water_training_rows':int((water_targets==0).sum()),
        'positive_water_training_share':float((water_targets>0).mean()),
        'mean_synthetic_water_fraction':float(water_targets.mean()),
        'mean_reference_water_fraction_diagnostic_only':float(refwater.mean()),
        'note':'The target-first positive-mixture prior is prescribed by the tutorial; reference prevalence is used for diagnosis only and never to alter training.'},indent=2))
    (out/'water_diagnostic_rules.json').write_text(json.dumps({'input_display_divisor':scale,'dark_quantile':q,
          'dark_VNIR_threshold':thresholds[0],'dark_SWIR_threshold':thresholds[1],
          'false_water_prediction_threshold':settings['water_false_prediction_threshold'],
          'water_absence_threshold':settings['water_absence_threshold'],
          'note':'Brightness is the mean of retained channels; darkness is a proxy, not an observed shadow label. No baseline correction.'},indent=2))
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    colours=['#d9654c','#27845b','#bb943a','#416a9d']
    for row,(sensor,pred) in enumerate(means.items()):
        signed=(pred[3,valid]-refwater)*100
        for col,(name,b,threshold) in enumerate([('VNIR',vnir[valid],thresholds[0]),('SWIR',swir[valid],thresholds[1])]):
            ax=axes[row,col]
            for target,class_name in enumerate(cfg['classes']):
                mask=dominant==target
                ax.scatter(b[mask],signed[mask],s=9,alpha=.45,color=colours[target],label=class_name)
            ax.axhline(0,color='black',lw=.7);ax.axvline(threshold,color='grey',lw=.8,ls=':')
            ax.set_title(f'{sensor} | {name} brightness');ax.set_xlabel('Mean input / 10,000');ax.set_ylabel('Predicted minus reference water (pp)')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Water overestimation diagnosis | reference used only to diagnose\nDotted line: image-only 20th percentile; dark is not a confirmed shadow')
    fig.savefig(out/'water_error_brightness.png',dpi=170,bbox_inches='tight');plt.close(fig)
    true=valid&(reference[3]>=.95)
    false=valid&(reference[3]<=settings['water_absence_threshold'])&(means['EnMAP'][3]>=settings['water_false_prediction_threshold'])
    median_brightness=np.median(brightness[valid_image])
    choices=[]
    if true.any():
        candidates=np.argwhere(true);values=brightness[true]
        choices.append(('median-brightness pure reference water',candidates[np.argmin(np.abs(values-np.median(values)))]))
    for label,mask in [('largest dark false-water prediction',false&shadow_like),
                        ('largest bright false-water prediction',false&(brightness>=median_brightness))]:
        if mask.any():
            candidates=np.argwhere(mask);choices.append((label,candidates[np.argmax(means['EnMAP'][3][mask])]))
    if len(choices)<2:
        candidates=np.argwhere(false)
        if len(candidates):choices.append(('largest false-water prediction',candidates[np.argmax(means['EnMAP'][3][false])]))
    with sqlite3.connect(root/cfg['data']['directory']/'library_berlin.gpkg') as con:
        names=[row[0] for row in con.execute('SELECT name FROM library_berlin ORDER BY fid')]
    library_water=np.flatnonzero(labels==3)
    fig,axes=plt.subplots(1,len(choices),figsize=(5*len(choices),4),layout='constrained',squeeze=False,sharey=True)
    examples=[]
    sections=np.split(np.arange(len(wavelengths)),np.flatnonzero(np.diff(wavelengths)>35)+1)
    for ax,(label,(r,c)) in zip(axes[0],choices):
        for part in sections:
            ax.plot(wavelengths[part],image[part,r,c]/scale,color='#b64d3a',lw=1.7,label='image pixel' if part[0]==0 else None)
            for j,index in enumerate(library_water):
                ax.plot(wavelengths[part],library[index,part]/scale,color=['#246a9c','#729dbb'][j%2],lw=1,ls=['--',':'][j%2],
                        label=names[index] if part[0]==0 else None)
        ax.set_title(label.replace(' ','\n',1),fontsize=10)
        ax.set_xlabel('Wavelength (nm)');ax.legend(fontsize=8)
        x,y=rasterio.transform.xy(profile['transform'],int(r),int(c))
        examples.append({'selection':label,'row':int(r),'col':int(c),'easting':x,'northing':y,
                         'reference_water':float(reference[3,r,c]),'predicted_water_enmap':float(means['EnMAP'][3,r,c]),
                         'dominant_reference':cfg['classes'][reference[:,r,c].argmax()],
                         'VNIR_brightness':float(vnir[r,c]),'SWIR_brightness':float(swir[r,c]),'shadow_like_proxy':bool(shadow_like[r,c])})
    axes[0,0].set_ylabel('Input / 10,000')
    fig.suptitle('Diagnostic examples versus both library water spectra | selected by fixed rules')
    fig.savefig(out/'water_example_spectra.png',dpi=170,bbox_inches='tight');plt.close(fig)
    pd.DataFrame(examples).to_csv(out/'water_example_pixels.csv',index=False)
    return groups,corr


def make_svr(settings,c,epsilon):
    return make_pipeline(StandardScaler(),SVR(C=c,epsilon=epsilon,gamma=settings['gamma'],kernel=settings['kernel'],cache_size=settings['cache_size_mb']))


def bounded(values,settings):
    return np.clip(values,0,1) if settings['clip_predictions_to_unit_interval'] else values


def select_svr_parameters(cfg,library,labels,weights,out):
    """Reference data cannot enter this interface or parameter-selection objective."""
    settings=cfg['stage_a_additions']['svr']
    validation_settings=dict(cfg['mixing'],mixtures_per_class=settings['validation_mixtures_per_class'],include_original_endmembers=False)
    datasets=[]
    for target,name in enumerate(cfg['classes']):
        seed=cfg['seed']+settings['selection_seed_offset']+target
        x,y=mixtures(library,labels,target,cfg['mixing'],seed)
        xv,yv=mixtures(library,labels,target,validation_settings,seed+1000)
        for sensor in ['EnMAP','S2-like']:
            datasets.append((name,sensor,x if sensor=='EnMAP' else x@weights,y,xv if sensor=='EnMAP' else xv@weights,yv))
    scores=[]
    for c,epsilon in product(settings['c_grid'],settings['epsilon_grid']):
        for name,sensor,x,y,xv,yv in datasets:
            model=make_svr(settings,c,epsilon);model.fit(x,y)
            score=np.mean(np.abs(bounded(model.predict(xv),settings)-yv))
            scores.append({'C':c,'epsilon':epsilon,'sensor':sensor,'class':name,'synthetic_validation_mae':score})
        LOG.info('SVR synthetic selection C=%s epsilon=%s pooled MAE=%.4f',c,epsilon,np.mean([s['synthetic_validation_mae'] for s in scores[-8:]]))
    scores=pd.DataFrame(scores);scores.to_csv(out/'svr_synthetic_parameter_scores.csv',index=False)
    ranking=scores.groupby(['C','epsilon'],as_index=False).synthetic_validation_mae.mean().sort_values(['synthetic_validation_mae','C','epsilon'])
    chosen=ranking.iloc[0].to_dict()
    chosen['selection']='Pooled synthetic validation MAE over four classes and both sensors; no reference access.'
    (out/'svr_selected_parameters.json').write_text(json.dumps(chosen,indent=2))
    return float(chosen['C']),float(chosen['epsilon'])


def svr_sensitivity(cfg,out,results,image,library,labels,weights,valid_image,reference,valid,profile):
    settings=cfg['stage_a_additions']['svr']
    c,epsilon=select_svr_parameters(cfg,library,labels,weights,out)
    n=cfg['mixing']['ensemble_size']
    pixels=image[:,valid_image].T
    sensor_x={'EnMAP':pixels,'S2-like':pixels@weights}
    members={sensor:np.zeros((n,4,len(pixels)),dtype='float32') for sensor in sensor_x}
    timings=[]
    for run in range(n):
        for target,name in enumerate(cfg['classes']):
            seed=cfg['seed']+run*100+target
            x,y=mixtures(library,labels,target,cfg['mixing'],seed)
            for sensor,values in sensor_x.items():
                tic=time.perf_counter();model=make_svr(settings,c,epsilon)
                model.fit(x if sensor=='EnMAP' else x@weights,y)
                chunk=settings['prediction_chunk_pixels'];outside=0
                for first in range(0,len(values),chunk):
                    raw=model.predict(values[first:first+chunk]);outside+=int(np.sum((raw<0)|(raw>1)))
                    members[sensor][run,target,first:first+chunk]=bounded(raw,settings)
                timing={'sensor':sensor,'run':run+1,'class':name,'seconds':time.perf_counter()-tic,
                        'raw_out_of_bounds_share':outside/len(values),'support_vectors':int(model[-1].n_support_.sum()),'seed':seed}
                timings.append(timing);LOG.info('SVR %s run %d %s: %.1fs (%d support vectors)',sensor,run+1,name,timing['seconds'],timing['support_vectors'])
    np.savez_compressed(results/'cache/svr_predictions.npz',**members)
    pd.DataFrame(timings).to_csv(out/'svr_runtime_and_bounds.csv',index=False)
    means,stds=maps_from_members(members,valid_image,valid.shape,n)
    summary,blocks,contexts=summarize(cfg,reference,valid,means,stds,n)
    diff,ci=differences(cfg,reference,valid,means,blocks,n)
    for name,table in [('svr_summary',summary),('svr_blocks',blocks),('svr_contexts',contexts),('svr_block_differences',diff),('svr_paired_bootstrap',ci)]:
        table.to_csv(out/f'{name}.csv',index=False)
    for sensor in means:
        stem=sensor.lower().replace('-','_')
        write_raster(out/f'svr_{stem}_fractions.tif',means[sensor],profile,cfg['classes'])
        write_raster(out/f'svr_{stem}_uncertainty.tif',stds[sensor],profile,cfg['classes'])
    return summary,ci


def run_additions(cfg,root,results):
    out=results/'stage_a_additions';out.mkdir(exist_ok=True)
    start=time.perf_counter()
    frozen=[results/'summary.csv',results/'cache/ensemble_predictions.npz']+list(results.glob('*_fractions_n*.tif'))
    before={path.name:sha256(path) for path in frozen}
    import yaml
    (out/'config_used.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False),encoding='utf-8')
    with np.load(results/'cache/prepared.npz') as data:
        image,library,labels=data['image'],data['library'],data['labels']
        wavelengths,weights=data['wavelengths'],data['weights']
        valid_image,reference,valid=data['valid_image'],data['reference'],data['valid']
    with rasterio.open(root/cfg['data']['directory']/'enmap_berlin.tif') as src:
        profile=src.profile.copy()
    with np.load(results/'cache/ensemble_predictions.npz') as data:
        members={sensor:data[sensor] for sensor in ['EnMAP','S2-like']}
    means,_=maps_from_members(members,valid_image,valid.shape,cfg['mixing']['ensemble_size'])
    ci,context_gain=gain_diagnostics(cfg,out,results,reference,valid,means)
    water_groups,water_corr=water_diagnosis(cfg,out,root,image,library,labels,wavelengths,valid_image,reference,valid,means,profile)
    svr,svr_ci=svr_sensitivity(cfg,out,results,image,library,labels,weights,valid_image,reference,valid,profile)
    rf=pd.read_csv(results/'summary.csv');rf=rf[rf.ensemble_size==cfg['mixing']['ensemble_size']].copy();rf['method']='RF'
    svr['method']='SVR'
    compare=pd.concat([rf,svr],ignore_index=True)
    compare.to_csv(out/'rf_svr_comparison.csv',index=False)
    write_additions_report(cfg,out,ci,context_gain,water_groups,water_corr,compare,svr_ci)
    after={path.name:sha256(path) for path in frozen}
    if before != after:
        raise RuntimeError('The frozen RF baseline changed during diagnostic additions.')
    (out/'rf_baseline_unchanged.json').write_text(json.dumps({'unchanged':True,'sha256':after},indent=2))
    (out/'runtime.json').write_text(json.dumps({'seconds':time.perf_counter()-start},indent=2))
    from report import write_report
    write_report(cfg,root,results,json.loads((results/'input_manifest.json').read_text(encoding='utf-8')),
                 pd.read_csv(results/'s2_response_audit.csv'),json.loads((results/'runtime.json').read_text(encoding='utf-8'))['stage_a_seconds'])
    LOG.info('Stage A additions complete in %.1fs',time.perf_counter()-start)


def write_additions_report(cfg,out,ci,context_gain,water_groups,water_corr,compare,svr_ci):
    gains=[]
    for method,frame in compare.groupby('method'):
        table=frame.pivot(index='class',columns='sensor',values='mae')
        table['gain']=table['S2-like']-table['EnMAP'];table['method']=method
        gains.append(table.reset_index())
    gains=pd.concat(gains,ignore_index=True)
    gains.to_csv(out/'rf_svr_sensor_gains.csv',index=False)
    rf_macro=gains[(gains.method=='RF')&(gains['class']=='macro')].iloc[0]
    svr_macro=gains[(gains.method=='SVR')&(gains['class']=='macro')].iloc[0]
    conclusion='the macro gain has the same sign' if rf_macro.gain*svr_macro.gain>0 else 'the macro conclusion changes sign'
    pivot=gains.pivot(index='class',columns='method',values='gain')
    sign_changes=pivot[(pivot.RF*pivot.SVR)<0].index.tolist()
    water=water_groups[(water_groups.sensor=='EnMAP')&(water_groups.context_type=='dominant_reference')]
    dark=water_groups[(water_groups.sensor=='EnMAP')&(water_groups.context_type=='shadow_like_proxy')]
    chosen=json.loads((out/'svr_selected_parameters.json').read_text(encoding='utf-8'))
    training_prior=json.loads((out/'water_synthetic_target_distribution.json').read_text(encoding='utf-8'))
    dark_rows=dark.set_index('context')
    text=f'''# Stage A additions — frozen RF baseline

## Does EnMAP improve MAE in this scene?

The existing **{cfg['evaluation']['bootstrap_resamples']} paired resamples of ~1 km blocks** give the 95% percentile intervals below. Positive differences mean lower EnMAP MAE. All 12 blocks with reference pixels enter the resampling with pixel-weighted aggregation; nine >=20-pixel blocks enter the supported-block share. Shares are unweighted proportions of blocks. This is within-scene evidence; neighbouring blocks, a shared parent image/library, and a small number of spatial units limit inference.

{markdown_table(ci[['class','s2_minus_enmap_mae','ci025','ci975','share_all_blocks_enmap_better','share_supported_blocks_enmap_better']])}

If an interval crosses zero, these data do not clearly resolve that class's gain under this resampling. A positive interval is not evidence for other scenes or missions. `rf_gain_by_context.csv` includes dominant-reference contexts, impervious gradients, and vegetation-mixture bins; the heatmap separates impervious/soil/vegetation confusions. Rare soil/water contexts have few pixels, and reference polygons are not spatially representative of the whole scene.

![Context gain](gain_by_context.png)

## Water overestimation: diagnosis only

The brightness thresholds, false-water rule and example-selection rules were fixed in `STAGE_A_ADDITIONS_DESIGN.md` before running these additions. Darkness requires low mean VNIR **and** SWIR relative to image-only 20th-percentile thresholds. It can represent water, shadow, dark materials or mixtures; no independent shadow labels exist here. The RF baseline remains unchanged, with no reference-trained correction and no sum-to-one experiment.

Water error by dominant reference context (EnMAP):

{markdown_table(water[['context','n_pixels','mean_predicted_water','mean_reference_water','water_bias','water_mae','false_water_pixel_share']])}

Darkness proxy:

{markdown_table(dark[['context','n_pixels','water_bias','water_mae','false_water_pixel_share']])}

Brightness associations (positive signed-error correlation means overprediction increases with brightness):

{markdown_table(water_corr[['sensor','brightness','subset','spearman_signed_water_error','spearman_absolute_water_error']])}

These associations diagnose confounding, not a causal effect of shadows or brightness. `water_example_pixels.csv` records deterministic true-water and false-water examples and their coordinates, fractions and brightness. The example spectra are diagnostic extrema/representative controls, not a random accuracy sample. Similar dark signatures and only two water endmembers limit spectral discrimination; systematic overestimation with low reference water is consistent with the tutorial's own diagnosis. The plots show whether the chosen false-water examples are dark or remain spectrally different from water; they do not by themselves establish the mechanism.

For this scene, water overprediction increases as brightness decreases (see the negative nonwater-only rank correlations), but it is **not confined to the darkness proxy**: {int(dark_rows.loc['other','false_water_pixels'])} false-water pixels are outside it versus {int(dark_rows.loc['dark_VNIR_and_SWIR','false_water_pixels'])} inside. The fixed training design is another plausible contributor: {training_prior['positive_water_training_share']:.1%} of water-target training rows contain positive water, with mean synthetic water fraction {training_prior['mean_synthetic_water_fraction']:.3f}, compared with mean reference water {training_prior['mean_reference_water_fraction_diagnostic_only']:.3f}. Only {training_prior['zero_water_training_rows']} original nonwater endmembers explicitly anchor zero water. This mismatch and the dark-surface spectral ambiguity are diagnostic explanations, not causal tests; no training prior was adjusted to the reference.

![Water brightness errors](water_error_brightness.png)
![Example spectra](water_example_spectra.png)

## SVR method sensitivity

The original methodological work used SVR ([Okujeni et al., 2017](https://doi.org/10.1109/JSTARS.2016.2634859)); the tutorial selects RF for speed. Our variant fits StandardScaler on **synthetic training mixtures only**, then uses RBF SVR. A single common **C={chosen['C']:g}, epsilon={chosen['epsilon']:g}, gamma=scale** was selected by pooled synthetic-validation MAE across all four targets and both sensors. The six candidate pairs and all synthetic scores are saved in `svr_synthetic_parameter_scores.csv`. Each target's 400 validation mixtures use separate seeds and exclude appended original endmembers; validation still shares the same material library and therefore measures mixture interpolation.

After selection, the pair was frozen and the three RF baseline seeds/mixtures were reused for paired SVR fits. Predictions are clipped to [0,1] by the predeclared sensitivity rule; raw out-of-range rates are in `svr_runtime_and_bounds.csv`. No sum-to-one constraint is applied, and no reference errors chose hyperparameters or bounds. This is a labelled SVR sensitivity, not a replacement RF baseline or a claim to replicate the original paper's full SVR setup.

{markdown_table(gains[['method','class','EnMAP','S2-like','gain']])}

For RF vs SVR, **{conclusion}**: macro S2-like-minus-EnMAP MAE is {rf_macro.gain:.5f} for RF and {svr_macro.gain:.5f} for SVR. Per-class signs and magnitudes must be inspected individually; estimator sensitivity does not validate a universal hyperspectral advantage.

Classes whose sensor-gain sign changes between RF and SVR: **{', '.join(sign_changes) if sign_changes else 'none'}**. This directly limits class-specific conclusions from the RF baseline.

SVR paired spatial-block intervals:

{markdown_table(svr_ci[['class','s2_minus_enmap_mae','ci025','ci975']])}

Full per-class accuracy, contexts, block metrics and uncertainty/error correlations are in `svr_summary.csv`, `svr_contexts.csv`, `svr_blocks.csv` and `rf_svr_comparison.csv`. The fixed RF baseline files remain unchanged.
'''
    (out/'REPORT.md').write_text(text,encoding='utf-8')
