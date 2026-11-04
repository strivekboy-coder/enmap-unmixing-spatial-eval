"""Publish measured Stage B trade-offs and refresh the project entry report."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import yaml

from report import markdown_table,write_report


def finalize(cfg,root,results):
    out=results/'stage_b'
    state=json.loads((out/'gate_status.json').read_text(encoding='utf-8'))
    if state['status']!='complete':return
    evidence=state['evidence']
    n_original=evidence['reference_support']['original_reference_pixels']
    n_common=evidence['reference_support']['common_supported_pixels']
    n_screen=evidence['reference_support']['consistency_screen_supported_pixels']
    n_excluded=evidence['reference_support']['cloud_feature_exclusions']
    year=evidence['selected_s2_year']
    days=evidence['clear_observation_stats']['days']
    summary=pd.read_csv(out/'summary.csv')
    full=summary[summary.evaluation_subset=='full_reference_supported']
    table=full.pivot(index='class',columns='sensor',values='mae').reindex(cfg['classes']+['macro'])
    counts=pd.DataFrame([{'class':name,'WorldCover_year':cfg['extension']['worldcover_year'],
                          'endmembers_used':evidence['endmembers_after_feature_filter'][name]} for name in cfg['classes']])
    confusion=pd.read_csv(out/'worldcover_reference_confusion.csv',index_col=0)
    class_agreement=[]
    for i,name in enumerate(cfg['classes']):
        total=int(confusion.loc[i].sum())
        agree=int(confusion.loc[i,str(i)]) if str(i) in confusion else 0
        class_agreement.append({'dominant_reference':name,'n_reference':total,'n_matching_WorldCover':agree,
                                 'agreement':agree/total if total else np.nan})
    agreement=pd.DataFrame(class_agreement)
    agreement.to_csv(out/'worldcover_agreement_by_reference_class.csv',index=False)
    c0=table['c0_S2_best_date'];c=table['c_S2_temporal'];d=table['d_late_fusion'];a=table['a_EnMAP_library']
    gain=pd.DataFrame({'class':table.index,'time_gain_c0_minus_c':(c0-c).values,
                       'fusion_gain_c_minus_d':(c-d).values,'fusion_gain_a_minus_d':(a-d).values})
    gain.to_csv(out/'model_gain_summary.csv',index=False)
    corr=full[['sensor','class','spearman_pixel_std_abs_error','spearman_block_std_mae','n_blocks']]
    correlation_table=corr.pivot(index='sensor',columns='class',values='spearman_pixel_std_abs_error')
    block_corr=corr[corr['class']=='macro'].set_index('sensor').spearman_block_std_mae
    correlation_table['macro_block_rho']=block_corr
    fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    with rasterio.open(out/'evaluation_common_support_mask.tif') as src:mask=src.read(1,masked=True).filled(0)>0
    with rasterio.open(results/'reference_fractions.tif') as src:reference=src.read(masked=True).filled(np.nan)
    for ax,model in zip(axes.ravel(),['a_EnMAP_library','b_S2_like_library','c0_S2_best_date','c_S2_temporal','d_late_fusion']):
        with rasterio.open(out/f'{model}_fractions.tif') as src:pred=src.read(masked=True).filled(np.nan)
        with rasterio.open(out/f'{model}_uncertainty.tif') as src:std=src.read(masked=True).filled(np.nan)
        err=np.abs(pred[:,mask]-reference[:,mask]).mean(axis=0)
        unc=std[:,mask].mean(axis=0)
        ax.hexbin(unc*100,err*100,gridsize=24,mincnt=1,bins='log',cmap='Blues')
        r=full[(full.sensor==model)&(full['class']=='macro')].iloc[0]
        ax.set_title(f'{model}\npixel rho={r.spearman_pixel_std_abs_error:.2f}; block rho={r.spearman_block_std_mae:.2f}',fontsize=10)
        ax.set_xlabel('Mean ensemble std across classes (pp)');ax.set_ylabel('Mean absolute fraction error (pp)')
    axes.ravel()[-1].axis('off')
    axes.ravel()[-1].text(.05,.8,f'Same {n_common:,} reference pixels\n{cfg["mixing"]["ensemble_size"]} ensemble members\n{int(full.n_blocks.iloc[0])} supported spatial blocks\n\nSpread is not calibrated confidence.\nKnown 2009 / {year} / {cfg["extension"]["worldcover_year"]} date gaps.',transform=axes.ravel()[-1].transAxes,va='top')
    fig.suptitle('Stage B uncertainty versus error | cross-year method demonstration')
    fig.savefig(out/'uncertainty_vs_error.png',dpi=170,bbox_inches='tight');plt.close(fig)
    attempts=json.loads((out/'seasonal_acquisition_proxy.json').read_text(encoding='utf-8'))
    note=f'''# Stage B results and limits

**This is a 2009 simulated-EnMAP / {year} real-S2 / {cfg['extension']['worldcover_year']} WorldCover method demonstration. The {year-2009}-year image gap and {cfg['extension']['worldcover_year']-year}-year S2/label gap prevent same-date sensor claims.**

All predefined gates passed. Both temporal models use the same external endmember pixels and the same 13-feature schema (one-date versus full-season statistics), with no reference-fraction training. The independent endmembers are outside the whole tutorial extent plus a 1 km gap.

{markdown_table(counts)}

## Accuracy on identical support

All five models use {n_common:,} reference pixels; the full original reference has {n_original:,}. {n_excluded} pixels lack common best-date/seasonal support. Under the fixed WorldCover consistency screen, {n_screen:,} common pixels remain. Units below are 0–1 fractions.

{markdown_table(table.reset_index())}

{markdown_table(gain)}

- Time-series macro gain over the same-source best date is **{(c0['macro']-c['macro'])*100:.3f} pp**. Water improves, while vegetation and impervious fractions worsen on this primary support. Thus the added temporal features do not uniformly close a spectral gap.
- Fixed fusion reduces macro MAE from **{c['macro']:.4f} to {d['macro']:.4f}** relative to temporal S2, and from **{a['macro']:.4f} to {d['macro']:.4f}** relative to the old EnMAP model.
- Fusion worsens water relative to temporal S2: **{c['water']:.4f} -> {d['water']:.4f}**. It also worsens soil relative to EnMAP alone: **{a['soil']:.4f} -> {d['soil']:.4f}**. Fixed averaging propagates each source's systematic errors; no weights were fitted to the reference to hide that trade-off.
- The "value of an acquisition" map is retrospective `MAE(c)-MAE(d)` against historical reference. Adding an old-date map is confounded by reference date and does not measure the benefit of obtaining a new EnMAP observation today.

## Consistency is not demonstrated stability

{markdown_table(agreement)}

**All 39 soil-dominated reference pixels fail the WorldCover agreement screen.** They remain in the primary full-support comparison, but none remains in the screened subset. The overall 81.4% agreement gate is dominated by vegetation/impervious locations and must not be read as class-wise agreement. WorldCover bare/sparse labels, historical soil definitions and real land-cover change cannot be separated here. The screened soil MAE is therefore not a test on soil-dominated sites. Lower screened error alone does not establish a more reliable model.

## Ensemble spread versus error

Pixel rank correlations by fraction and macro block correlation (only nine supported blocks):

{markdown_table(correlation_table.reset_index())}

These are associations under model/synthetic-mixture variability, not calibrated uncertainty or independent-scene evidence.

![Uncertainty versus error](uncertainty_vs_error.png)

## Acquisition proxy and data audit

The strict local usable-area rule succeeds on {int(round(attempts['usable_probability']*days))}/{days} dates in {year}: plug-in probability **{attempts['usable_probability']:.4f}**, giving **{attempts['expected_attempts']:.2f}** attempts under an independent, stationary S2-to-EnMAP cloudiness proxy. SCL classes 4/5/6 also reflect classification/shadow behaviour, and zero-success months do not imply a truly zero probability. The separate scene-level <=20% cloud proxy is in `monthly_acquisition_proxy.csv`; it is not interchangeable with local SCL usability.

Five selected assets had contradictory offset metadata. Same-date red/NIR/SWIR comparisons with legacy assets, without using reference labels or prediction errors, establish that their -0.1 correction had already been applied. The explicit item whitelist and audit guard prevent double subtraction. Declared and effective offsets, dates, item IDs and URLs are saved in `radiometric_calibration_by_date_band.csv`. This is a data-quality correction, not a fitted model adjustment.

Full methods and both support tables: [STATUS.md](STATUS.md). Core evidence: `summary.csv`, `blocks.csv`, `contexts.csv`, `model_gain_summary.csv`, `value_of_acquisition_blocks.csv`, `evaluation_support_by_context.csv`, and the two mask GeoTIFFs. Summary-feature mixing of NDVI/medians/amplitude is an approximation rather than physically exact mixing of complete time series.

![Retrospective acquisition value](value_of_acquisition.png)

Stage C scaling, official QGIS comparison and real EnMAP transfer remain later additions.
'''
    (out/'RESULTS_SUMMARY.md').write_text(note,encoding='utf-8')
    status_text=(out/'STATUS.md').read_text(encoding='utf-8')
    if 'External endmembers retained after fixed cloud/feature filtering:' not in status_text:
        status_text=status_text.replace('## Primary matched-support results',
            'External endmembers retained after fixed cloud/feature filtering:\n\n'+markdown_table(counts)+'\n\n## Primary matched-support results')
    if '## Reference consistency by class' not in status_text:
        status_text=status_text.replace('## Method and provenance',
            '## Reference consistency by class\n\n'+markdown_table(agreement)+'\n\n**No soil-dominated reference pixels remain in the consistency-screened subset.** Use the full-support comparison as primary; the screen does not establish temporal stability. See [measured trade-offs](RESULTS_SUMMARY.md).\n\n## Method and provenance')
    (out/'STATUS.md').write_text(status_text,encoding='utf-8')
    write_report(cfg,root,results,json.loads((results/'input_manifest.json').read_text(encoding='utf-8')),
                 pd.read_csv(results/'s2_response_audit.csv'),json.loads((results/'runtime.json').read_text(encoding='utf-8'))['stage_a_seconds'])
    print('Stage B reports/uncertainty figure refreshed; class-wise zero soil-screen support is explicit.')


if __name__=='__main__':
    root=Path(__file__).resolve().parent
    cfg=yaml.safe_load((root/'config.yaml').read_text(encoding='utf-8'))
    finalize(cfg,root,root/cfg['results'])
