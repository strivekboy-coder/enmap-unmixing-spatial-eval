"""Inspect the provided data before choosing any preprocessing parameters."""
import json
import sqlite3
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio

root = Path(__file__).parent
with rasterio.open(root / 'data/enmap_berlin.tif') as src:
    print('PROFILE', src.profile)
    print('TAGS', src.tags())
    print('BAND 1', src.tags(1), src.descriptions[:3])
    print('ENVI', src.tags(ns='ENVI'))
    image = src.read()
    print('RANGE', np.nanmin(image), np.nanmax(image), 'ZERO', (image == 0).mean())
with sqlite3.connect(root / 'data/library_berlin.gpkg') as con:
    profiles = [json.loads(row[0]) for row in con.execute('SELECT profiles FROM library_berlin')]
    print('LIBRARY KEYS', profiles[0].keys())
    print('LIBRARY AXIS', {k: v for k, v in profiles[0].items() if k != 'y'})
    print('LIBRARY Y', profiles[0]['y'][:5])
    print('CLASSES', con.execute('SELECT level_1,count(*) FROM library_berlin GROUP BY level_1').fetchall())
ref = gpd.read_file(root / 'data/landcover_berlin.gpkg')
print('REFERENCE', ref.crs, ref.columns.tolist(), ref.groupby('level_1').size().to_dict(), ref.total_bounds)
