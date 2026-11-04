# Real EnMAP first look — Berlin, 2025-05-11

**16-year calendar gap:** 2025 real EnMAP versus the 2009 simulated scene/reference.
Historical discrepancies combine model/sensor, geometry, season and genuine
land-cover change; they do not establish same-date accuracy.

XML confirms centre acquisition **2025-05-11 10:42:03.869622 UTC**, processing-chain
version **01.05.02**, 224 bands (418.416–2445.3 nm),
viewing zenith **9.557 degrees**,
viewing azimuth **103.265 degrees**,
solar elevation/azimuth **55.329/171.341 degrees**.
Reflectance is **DN x 0.0001 + 0**, applied exactly once; raw nodata is **-32768**.
Land mode, bilinear, terrain correction, no cirrus/haze removal or band interpolation.

Atmospheric/testflag-clear nearest-grid share: 86.60%. **Conservative complete-feature inference support: 2,721/88,000 pixels (3.09% of the tutorial rectangle).** Historical reference overlap: 31 pixels. Reject cloud, cloud shadow, haze, any cirrus, snow, defective source bands and any quality-test flag; every bilinear contributor must pass. Individual flag shares overlap.

Historical >=95%-pure atmospheric-clear support: roof 11, low vegetation 19, soil 14, water 8. Labels describe 2009, not independently verified 2025 materials. Water-vapour windows and roughness/finite-support flags (not measured SNR) are marked. The historical soil and water groups have 2025 median-spectrum NDVI 0.50 and 0.60; so these are not verified contemporary pure endmembers.

Frozen 3-member library RF applied to 2,721 real QA-clear pixels. MAE below uses the identical 31 reference pixels for both dates.

| Fraction | Real 2025 vs 2009 MAE | Simulated 2009 MAE, same pixels |
| --- | ---: | ---: |
| impervious | 0.1305 | 0.1046 |
| vegetation | 0.2347 | 0.0946 |
| soil | 0.1647 | 0.0884 |
| water | 0.2081 | 0.2114 |
| macro | 0.1845 | 0.1248 |

Dominant-reference support (impervious/vegetation/soil/water): 17/10/1/3; these 31 spatially sparse pixels cannot establish generalization.

Alignment: bilinear spatial resampling; linear 177-band interpolation without
extrapolation/gap bridging. One edge maps 1783 to 1780.22 nm (within half-FWHM/5nm).
No reference-fitted correction, tuning or sum-to-one. All 12 rebuilt EnMAP models
reproduce Stage A predictions exactly. S2-like is sensitivity-only: 2/71,158 replay
pixels differ, maximum 0.009905. Three-member SD is not calibrated uncertainty.

**Actual-S2 check:** Sentinel-2C 2025-05-12 (+1 day), scene cloud 8.00%. AROSICS 1.14.0 accepted 1 independent-mask 64x64 window (reliability 68.95%, an internal score). Correction to add to EnMAP origin: east +4.85 m, north -7.39 m. This is a local estimate, not scene-wide accuracy; RF results remain unshifted. Original ESA XML and red/NIR/SWIR JP2 comparisons (59,755 pixels/band) establish that S2 COG DN already includes -1000 DN: effective reflectance offset is 0, avoiding double subtraction. No EnMAP-based radiometric fitting. Official S2C SRFs give 10 band comparisons, 10,753–47,727 jointly qualified pairs/band: bias (EnMAP minus S2) -0.0267 to -0.0090; RMSD 0.0259–0.0514. One-day, atmosphere, view/sun, PSF and resampling differences prevent calibration-validation claims.

Raw data/models/tables/GeoTIFFs stay in ignored `data/real_enmap/20250511/`;
credentials never enter Git. Deadline: **2026-11-04 15:30 Europe/Berlin**.
Re-run Stage A first, then `python real_enmap.py --stage STAGE` in this order:
metadata/download/prepare-models/quality/spectra/transfer; optional-prepare/optional-audit/optional-coreg/optional-compare.
Use `real_enmap_config.yaml`; set a new explicit deadline for later reruns.

[DLR scene metadata](https://geoservice.dlr.de/eoc/ogc/stac/v1/collections/ENMAP_HSI_L2A/items/ENMAP01-____L2A-DT0000129412_20250511T104201Z_009_V010502_20250515T040235Z) ·
[Quality codes: product specification section 4.5](https://www.enmap.org/data/doc/EN-PCV-ICD-2009-2_HSI_Product_Specification_Level1_Level2.pdf).

[spectra_quality.png](results/real_data_first_look/spectra_quality.png)

[transfer_summary.png](results/real_data_first_look/transfer_summary.png)
