"""Verify saved scientific artifacts against independently read GeoTIFF values."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml

from report import write_report

root = Path(__file__).parent
cfg = yaml.safe_load((root/'config.yaml').read_text(encoding='utf-8'))
out = root/cfg['results']
summary = pd.read_csv(out/'summary.csv')
with rasterio.open(out/'reference_fractions.tif') as src:
    ref = src.read(masked=True)
    valid = ~np.ma.getmaskarray(ref).any(axis=0)
    observed = ref.data[:, valid]
assert valid.sum() == 1578
for size in [3, 5]:
    for sensor, stem in [('EnMAP', 'enmap'), ('S2-like', 's2_like')]:
        with rasterio.open(out/f'{stem}_fractions_n{size}.tif') as src:
            pred = src.read(masked=True)[:, valid]
            assert src.descriptions == tuple(cfg['classes'])
            assert src.crs.to_epsg() == 32633
            assert np.isfinite(pred).all()
            assert np.all((pred >= 0) & (pred <= 1))
        with rasterio.open(out/f'{stem}_uncertainty_n{size}.tif') as src:
            std = src.read(masked=True)[:, valid]
            assert np.all(std >= 0)
        table = summary[(summary.sensor == sensor) & (summary.ensemble_size == size)].set_index('class')
        for i, name in enumerate(cfg['classes']):
            error = np.asarray(pred[i])-observed[i]
            np.testing.assert_allclose(np.abs(error).mean(), table.loc[name, 'mae'], atol=1e-7)
            np.testing.assert_allclose(np.sqrt((error**2).mean()), table.loc[name, 'rmse'], atol=1e-7)
for name in ['fraction_maps', 'block_mae_difference', 'uncertainty_vs_error', 'spectral_importance', 'acquisition_priority']:
    assert (out/f'{name}.png').stat().st_size > 10000
manifest = json.loads((out/'input_manifest.json').read_text(encoding='utf-8'))
srf = pd.read_csv(out/'s2_response_audit.csv')
runtime = json.loads((out/'runtime.json').read_text(encoding='utf-8'))
write_report(cfg, root, out, manifest, srf, runtime['stage_a_seconds'])
print('Verified all 40 saved class/sensor/ensemble metric combinations against GeoTIFFs; five figures present. README refreshed.')
