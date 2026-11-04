# Real EnMAP first look — Berlin, 2025-05-11

**Preliminary transfer experiment with a 16-year calendar gap:** real EnMAP (2025)
versus the 2009 simulated scene/reference. Later differences include genuine land
cover change, season, radiometry, geometry and model error; they cannot establish
same-date accuracy.

Scene: `DT0000129412`, tile `009`. XML was downloaded and parsed successfully.

| XML field | Measured metadata |
| --- | --- |
| Scene-centre acquisition | 2025-05-11 10:42:03.869622 UTC |
| Processing-chain / archived L0 version | 01.05.02 / 01.05.02 |
| Viewing zenith / azimuth | 9.556863 / 103.264674 degrees |
| Sun elevation / azimuth | 55.328669 / 171.340563 degrees |
| Spectral bands | 224; 418.416–2445.3 nm |
| Reflectance conversion | DN x 0.0001 + 0, every band, applied once |
| Raw invalid value | -32768 |
| Processing | Land mode; bilinear; no cirrus/haze removal; no band interpolation |

**Progress:** step 1 complete. Quality layers have arrived; image download and
steps 2–4 are underway. No tutorial-area clear fraction, real spectrum or transfer
accuracy is claimed yet. Optional AROSICS/Sentinel-2 checks have not started.

Stage A estimator objects were not originally saved. All 12 reconstructed EnMAP
estimators exactly reproduce the saved 2009 predictions. One reconstructed
S2-like estimator differs at 2/71,158 simulated pixels (maximum fraction difference
0.009905); that branch is marked sensitivity-only. No reference-based tuning occurs.

Raw downloads, XML, masks, model objects, numerical diagnostics and GeoTIFFs stay
in ignored `data/real_enmap/20250511/`. Only code/configuration, up to two summary
figures and this report are committed. Credentials are read locally, never logged.
Geoservice authentication now uses the IPS main identity, as specified by DLR.

Hard stop: **2026-11-04 15:30 Europe/Berlin**. Run
`python real_enmap.py --stage metadata` then `--stage download` with
`real_enmap_config.yaml`. Network work and model execution have deadline guards.

Sources: [scene metadata](https://geoservice.dlr.de/eoc/ogc/stac/v1/collections/ENMAP_HSI_L2A/items/ENMAP01-____L2A-DT0000129412_20250511T104201Z_009_V010502_20250515T040235Z),
[product specification, quality codes in section 4.5](https://www.enmap.org/data/doc/EN-PCV-ICD-2009-2_HSI_Product_Specification_Level1_Level2.pdf),
[DLR IPS identity guidance](https://www.enmap.org/news/2026-10-16_2/).
