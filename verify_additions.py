"""Audit SVR artifacts, parameter selection, block shares and frozen RF outputs."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml

from data_io import sha256
from stage_a_additions import maps_from_members,water_diagnosis,write_additions_report

root=Path(__file__).parent
cfg=yaml.safe_load((root/'config.yaml').read_text(encoding='utf-8'))
results=root/'results';out=results/'stage_a_additions'
frozen=json.loads((out/'rf_baseline_unchanged.json').read_text(encoding='utf-8'))
for name,digest in frozen['sha256'].items():
    path=results/'cache'/name if name=='ensemble_predictions.npz' else results/name
    assert sha256(path)==digest,f'RF artifact changed: {name}'
scores=pd.read_csv(out/'svr_synthetic_parameter_scores.csv')
rank=scores.groupby(['C','epsilon']).synthetic_validation_mae.mean().sort_values()
chosen=json.loads((out/'svr_selected_parameters.json').read_text(encoding='utf-8'))
assert tuple(rank.index[0])==(chosen['C'],chosen['epsilon'])
assert scores.groupby(['C','epsilon']).size().eq(8).all()
with np.load(results/'cache/prepared.npz') as data:
    image,library,labels,wavelengths=data['image'],data['library'],data['labels'],data['wavelengths']
    valid_image,reference,valid=data['valid_image'],data['reference'],data['valid']
with rasterio.open(root/'data/enmap_berlin.tif') as src:profile=src.profile.copy()
svr=pd.read_csv(out/'svr_summary.csv')
for sensor in ['EnMAP','S2-like']:
    stem=sensor.lower().replace('-','_')
    with rasterio.open(out/f'svr_{stem}_fractions.tif') as src:pred=src.read(masked=True)[:,valid]
    rows=svr[svr.sensor==sensor].set_index('class')
    for target,name in enumerate(cfg['classes']):
        np.testing.assert_allclose(np.abs(np.asarray(pred[target])-reference[target,valid]).mean(),rows.loc[name,'mae'],atol=1e-7)
    assert ((pred>=0)&(pred<=1)).all()
blocks=pd.read_csv(out/'rf_block_gains_all_support.csv')
ci=pd.read_csv(out/'rf_gain_intervals_and_block_shares.csv')
assert ci.bootstrap_resamples.min()>=1000
for _,row in ci.iterrows():
    group=blocks[blocks['class']==row['class']]
    np.testing.assert_allclose(group.enmap_better.mean(),row.share_all_blocks_enmap_better)
    np.testing.assert_allclose(group[group.passes_block_support].enmap_better.mean(),row.share_supported_blocks_enmap_better)
with np.load(results/'cache/ensemble_predictions.npz') as data:
    means,_=maps_from_members({s:data[s] for s in ['EnMAP','S2-like']},valid_image,valid.shape,3)
water,corr=water_diagnosis(cfg,out,root,image,library,labels,wavelengths,valid_image,reference,valid,means,profile)
write_additions_report(cfg,out,ci,pd.read_csv(out/'rf_gain_by_context.csv'),water,corr,
                        pd.read_csv(out/'rf_svr_comparison.csv'),pd.read_csv(out/'svr_paired_bootstrap.csv'))
print('Verified frozen RF checksums, synthetic-only parameter ranking, all eight SVR class metrics and paired block shares; refreshed diagnosis/report without refitting.')
