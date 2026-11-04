"""Compare reprocessed S2 assets with same-date legacy assets, without reference labels."""
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from temporal import choose_daily_items

logging.getLogger('rasterio.session').setLevel(logging.WARNING)
root=Path(__file__).parent
out=root/'results/stage_b'
features=json.loads((out/'cube/stac_2018.json').read_text(encoding='utf-8'))['features']
with rasterio.open(root/'data/enmap_berlin.tif') as src:
    profile=src.profile.copy()
    valid=src.read_masks(1)>0
from rasterio.warp import transform_bounds
bbox=transform_bounds(profile['crs'],'EPSG:4326',*rasterio.transform.array_bounds(profile['height'],profile['width'],profile['transform']))
selected=choose_daily_items(features,bbox)

def read_band(item,band):
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_TIMEOUT=60):
        with rasterio.open(item['assets'][band]['href']) as src:
            tags={'gdal_scale':src.scales[0],'gdal_offset':src.offsets[0]}
            with WarpedVRT(src,crs=profile['crs'],transform=profile['transform'],width=profile['width'],height=profile['height'],
                           resampling=Resampling.average,nodata=0,dtype='uint16') as vrt:
                return vrt.read(1).astype('float32'),tags

records=[]
for item in selected:
    meta=item['assets']['red']['raster:bands'][0]
    if meta.get('offset',0)==0:continue
    day=item['properties']['datetime'][:10]
    tile=item['id'].split('_')[1]
    candidates=[x for x in features if x['properties']['datetime'][:10]==day and x['id'].split('_')[1]==tile
                and x['assets']['red']['raster:bands'][0].get('offset',0)==0]
    if not candidates:
        records.append({'date':day,'selected_item':item['id'],'status':'no_same_date_legacy_match'})
        continue
    legacy=sorted(candidates,key=lambda x:x['id'])[0]
    for band in ['red','nir','swir16']:
        meta=item['assets'][band]['raster:bands'][0]
        a,at=read_band(item,band);b,bt=read_band(legacy,band)
        ok=valid&(a>0)&(b>0)
        oldmeta=legacy['assets'][band]['raster:bands'][0]
        delta_raw=float(np.median((a-b)[ok]))
        delta_cal=float(np.median((a*meta['scale']+meta['offset']-b*oldmeta['scale']-oldmeta['offset'])[ok]))
        records.append({'date':day,'band':band,'selected_item':item['id'],'legacy_item':legacy['id'],'n_pixels':int(ok.sum()),
                        'median_raw_DN_difference':delta_raw,'median_declared_calibrated_difference':delta_cal,
                        'selected_scale':meta['scale'],'selected_offset':meta['offset'],
                        'selected_flag':item['properties'].get('earthsearch:boa_offset_applied'),
                        'selected_negative_reflectance_share':float(((a*meta['scale']+meta['offset'])[ok]<0).mean()),
                        **at,'status':'same_date_comparison'})
        print(records[-1],flush=True)
pd.DataFrame(records).to_csv(out/'same_date_radiometric_audit.csv',index=False)
