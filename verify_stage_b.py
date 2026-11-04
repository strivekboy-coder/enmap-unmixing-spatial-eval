"""Verify external separation, common supports, calibration and saved fusion outputs."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml

root=Path(__file__).parent
cfg=yaml.safe_load((root/'config.yaml').read_text(encoding='utf-8'))
out=root/'results/stage_b'
status=json.loads((out/'gate_status.json').read_text(encoding='utf-8'))
if status['status']=='stopped':
    assert status.get('reason')
    print('Stage B stopped at a documented fixed gate:',status['reason'])
    raise SystemExit(0)
assert status['status']=='complete',status
evidence=status['evidence']
assert evidence['worldcover_reference_agreement']>=cfg['extension']['minimum_worldcover_reference_agreement']
assert min(evidence['endmembers_after_feature_filter'].values())>=30
years=pd.read_csv(out/'year_selection.csv').sort_values(['median_clear_dates','p10_clear_dates','year'],ascending=[False,False,True])
assert years.iloc[0].year==evidence['selected_s2_year']
assert years.iloc[0].supported_share>=cfg['extension']['minimum_supported_pixel_share']
endmembers=pd.read_csv(out/'external_endmembers_used.csv')
assert (endmembers.worldcover_purity>=cfg['extension']['minimum_endmember_purity']).all()
assert (endmembers.clear_observation_count>=5).all()
with rasterio.open(root/'data/enmap_berlin.tif') as src:bounds=src.bounds
dx=np.maximum(np.maximum(bounds.left-endmembers.x,endmembers.x-bounds.right),0)
dy=np.maximum(np.maximum(bounds.bottom-endmembers.y,endmembers.y-bounds.top),0)
distance=np.hypot(dx,dy)
assert (distance>cfg['extension']['buffer_inner_m']).all()
assert (distance<=cfg['extension']['buffer_outer_m']).all()
cal=pd.read_csv(out/'radiometric_calibration_by_date_band.csv')
corrected=cal[cal.audited_override]
assert set(corrected['item'])==set(cfg['extension']['audited_offset_already_applied_items'])
np.testing.assert_allclose(corrected.offset,0)
np.testing.assert_allclose(corrected.declared_offset,-.1)
summary=pd.read_csv(out/'summary.csv')
with rasterio.open(root/'results/reference_fractions.tif') as src:reference=src.read(masked=True).filled(np.nan)
masks={}
for name,file in [('full_reference_supported','evaluation_common_support_mask.tif'),
                   ('worldcover_consistent_supported','evaluation_consistency_support_mask.tif')]:
    with rasterio.open(out/file) as src:masks[name]=src.read(1,masked=True).filled(0)>0
    assert summary[summary.evaluation_subset==name].n_pixels.nunique()==1
    assert int(summary[summary.evaluation_subset==name].n_pixels.iloc[0])==int(masks[name].sum())
models={}
for model in summary.sensor.unique():
    with rasterio.open(out/f'{model}_fractions.tif') as src:pred=src.read(masked=True).filled(np.nan)
    models[model]=pred
    for subset,mask in masks.items():
        rows=summary[(summary.sensor==model)&(summary.evaluation_subset==subset)].set_index('class')
        for target,name in enumerate(cfg['classes']):
            error=pred[target,mask]-reference[target,mask]
            np.testing.assert_allclose(np.abs(error).mean(),rows.loc[name,'mae'],atol=1e-7)
            np.testing.assert_allclose(np.sqrt((error**2).mean()),rows.loc[name,'rmse'],atol=1e-7)
mask=masks['full_reference_supported']
np.testing.assert_allclose(models['d_late_fusion'][:,mask],
                           (models['a_EnMAP_library'][:,mask]+models['c_S2_temporal'][:,mask])/2,atol=1e-7)
with np.load(out/'cube/temporal_model_members.npz') as members:
    np.testing.assert_allclose(members['d_late_fusion'],(members['a_EnMAP_library']+members['c_S2_temporal'])/2,atol=1e-7,equal_nan=True)
assert evidence['s2_minus_2009_years']==9
features=pd.read_csv(out/'external_endmember_features.csv')
np.testing.assert_allclose(features.ndvi_one_date_median,features.ndvi_one_date_max)
np.testing.assert_allclose(features.ndvi_one_date_amplitude,0)
print('Verified independent external locations, all gates, two common evaluation masks, 40 saved class metrics, five audited offset corrections and fixed member-wise late fusion.')
