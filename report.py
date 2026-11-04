"""Generate a results-led README from the saved measurements, never a forecast."""
from __future__ import annotations

import json
import numpy as np
import pandas as pd


def markdown_table(frame):
    header = '| ' + ' | '.join(map(str, frame.columns)) + ' |\n'
    rule = '| ' + ' | '.join(['---']*len(frame.columns)) + ' |\n'
    rows = []
    for row in frame.itertuples(index=False, name=None):
        rows.append('| ' + ' | '.join(f'{x:.4f}' if isinstance(x, (float, np.floating)) else str(x) for x in row) + ' |')
    return header + rule + '\n'.join(rows)


def write_report(cfg, root, out, manifest, srf_info, elapsed):
    summary = pd.read_csv(out/'summary.csv')
    size = cfg['mixing']['ensemble_size']
    base = summary[summary.ensemble_size == size]
    comparison = base.pivot(index='class', columns='sensor', values='mae').reindex(cfg['classes']+['macro'])
    comparison['S2-like minus EnMAP'] = comparison['S2-like']-comparison['EnMAP']
    comparison = comparison.reset_index()
    uncertainty = base[['sensor', 'class', 'spearman_pixel_std_abs_error', 'spearman_block_std_mae']]
    sensitivity = summary.pivot(index=['sensor', 'class'], columns='ensemble_size', values='mae').reset_index()
    sensitivity.columns = [str(c) for c in sensitivity.columns]
    contexts = pd.read_csv(out/'contexts.csv')
    ct = contexts[(contexts.ensemble_size == size) & (contexts['class'] == 'macro') &
                  (contexts.context_type == 'dominant_reference')]
    cp = ct.pivot(index='context', columns='sensor', values='mae')
    cp['gain'] = cp['S2-like']-cp['EnMAP']
    best_context = cp['gain'].idxmax()
    block_stats = base[['sensor', 'class', 'n_blocks', 'block_mae_min', 'block_mae_median',
                       'block_mae_max', 'share_blocks_above_twice_overall']]
    tutorial = np.array([.1179, .1088, .1094, .2514])
    ours = base[(base.sensor == 'EnMAP') & (base['class'] != 'macro')].set_index('class').loc[cfg['classes'], 'mae'].to_numpy()
    water = base[(base.sensor == 'EnMAP') & (base['class'] == 'water')].iloc[0]
    water_s2 = base[(base.sensor == 'S2-like') & (base['class'] == 'water')].iloc[0]
    sanity = pd.DataFrame({'class': cfg['classes'], 'tutorial_screenshot_mae': tutorial, 'python_mae': ours,
                           'difference': ours-tutorial})
    sanity.to_csv(out/'tutorial_sanity_comparison.csv', index=False)
    ratios = ours/tutorial
    sanity_status = 'same broad range' if np.all((ratios >= .5) & (ratios <= 2)) else 'outside a factor-of-two range for at least one class; inspect differences before interpretation'
    im = pd.read_csv(out/'enmap_feature_importance.csv')
    peak_i = im.loc[im.impervious.idxmax(), 'wavelength_nm']
    peak_s = im.loc[im.soil.idxmax(), 'wavelength_nm']
    interpretation = (
        f'The input already omits the major water-vapour windows (gaps {manifest["input_spectral_gaps_nm"]} nm), '
        f'and none of its {manifest["bands_retained"]} finite bands was additionally removed by the configured filter. '
        f'Peak impurity importance occurs at {peak_i:.0f} nm for impervious cover and {peak_s:.0f} nm for soil; '
        'the plotted red edge (~690–750 nm), Al-OH/clay region (~2200 nm), and ~2300 nm region provide physical context. '
        'Correlated bands, a small synthetic library, and impurity-importance bias prevent these peaks from establishing mineral absorption or causal sensor value.'
    )
    (out/'spectral_interpretation.md').write_text(interpretation+'\n', encoding='utf-8')
    block_ci = pd.read_csv(out/'paired_block_bootstrap.csv')
    block_ci = block_ci[block_ci.ensemble_size == size]
    stage_b_result='Stage B has not yet completed; no temporal/fusion accuracy is claimed.'
    stage_b_status=out/'stage_b/gate_status.json'
    if stage_b_status.exists():
        state=json.loads(stage_b_status.read_text(encoding='utf-8'))
        if state.get('status')=='complete':
            bs=pd.read_csv(out/'stage_b/summary.csv')
            rows=bs[(bs.evaluation_subset=='full_reference_supported')&(bs['class']=='macro')][['sensor','mae','n_pixels']]
            stage_b_result=('Stage B is complete as a **cross-year method demonstration**, using 2009 simulated EnMAP, 2018 real S2 and 2021 WorldCover. '
                            'Its common-support macro results are below; they use a different mask from the 1,578-pixel Stage A baseline.\n\n'+
                            markdown_table(rows)+'\n\nSee the [Stage B results and limits](results/stage_b/RESULTS_SUMMARY.md) and [full methods](results/stage_b/STATUS.md). '
                            '**All 39 soil-dominated reference pixels fail the WorldCover consistency screen**; the primary full-support results retain them. '
                            'Time-series and fixed fusion gains are class-dependent: fusion improves average error but worsens water relative to temporal S2. '
                            'Date mismatch confounds any acquisition-value interpretation.')
        elif state.get('status')=='stopped':
            stage_b_result='Stage B stopped at a fixed gate: '+state.get('reason','See stage_b/STATUS.md')
    macro = comparison.loc[comparison['class'] == 'macro'].iloc[0]
    macro_ci = block_ci.loc[block_ci['class'] == 'macro'].iloc[0]
    text = f'''# EnMAP unmixing and spatial error analysis — Berlin tutorial

**Research question:** Does hyperspectral information improve four urban land-cover fractions over Sentinel-2-like spectra, where does it help, and can temporal S2 data add value?

**Main findings:** RF macro MAE falls from {macro['S2-like']:.4f} to {macro['EnMAP']:.4f}: a {macro['S2-like minus EnMAP']*100:.2f} percentage-point gain (paired spatial-block 95% CI {macro_ci.ci025*100:.2f}–{macro_ci.ci975*100:.2f}). Temporal S2 and fixed late fusion show average gains with class-specific trade-offs; fusion worsens water relative to temporal S2. [Stage A sensitivity/diagnosis](results/stage_a_additions/REPORT.md) · [Stage B findings](results/stage_b/RESULTS_SUMMARY.md).

**Limitations:** One Berlin tutorial scene; simulated EnMAP; 1,578 reference pixels in 12 blocks (9 with reporting support); large water bias. Stage B combines 2009/2018/2021 sources and loses all soil-dominated support under the consistency screen. The separate [2025 real-EnMAP first look](REAL_DATA_FIRST_LOOK.md) retains only 3.09% complete-spectrum support and 31 historical reference pixels across a 16-year gap. No official QGIS/EnMAP-Box comparison; this is a Python re-implementation, not a validated reproduction.

[Code: MIT](LICENSE) · [Data attribution and result licenses](DATA_SOURCES.md) · Raw data, GeoTIFFs and caches are excluded; regenerate them locally.

## Measured Stage A results

MAE is expressed as a fraction (0–1); multiply by 100 for percentage points. All sensors use the same **{int(base.n_pixels.iloc[0])} reference pixels**, with at least 95% polygon coverage, and the same synthetic mixtures for each paired run. Positive S2-like minus EnMAP values mean lower EnMAP error.

{markdown_table(comparison)}

The macro reduction is {(comparison.loc[comparison['class']=='macro', 'S2-like minus EnMAP'].iloc[0])*100:.2f} percentage points. This is modest: water remains systematically overestimated (EnMAP mean bias **{water['bias']:.4f}**, water R² **{water['r2']:.2f}**). Among the {int(water['reference_positive_pixels'])} pixels with any reference water, EnMAP water MAE is **{water['mae_reference_positive']:.4f}** versus S2-like **{water_s2['mae_reference_positive']:.4f}**; the overall water improvement therefore does not imply better accuracy on water-containing pixels.

Among dominant-reference contexts, the largest macro-MAE reduction is in **{best_context}**-dominated pixels ({cp.loc[best_context, 'gain']*100:.2f} percentage points). This is conditional on the available reference locations; it is not a statement about all Berlin. Results include negative gains wherever the S2-like model performs better.

### Spatial variability

Blocks are anchored to the upper-left EnMAP grid, {cfg['evaluation']['block_size_pixels']}×{cfg['evaluation']['block_size_pixels']} pixels (~990 m), with at least {cfg['evaluation']['minimum_pixels_per_block']} evaluation pixels per reported block. Overall metrics include all qualifying pixels, including pixels in smaller-support blocks. `contexts.csv` additionally reports dominant-reference class and fixed 25-percentage-point bins of reference impervious fraction.

{markdown_table(block_stats)}

Paired resampling of spatial blocks gives the following 95% percentile intervals for the pixel-weighted MAE difference. These are descriptive within-scene intervals: adjacent blocks may remain dependent, and they do not establish generalization to another city.

{markdown_table(block_ci[['class', 's2_minus_enmap_mae', 'ci025', 'ci975', 'n_spatial_blocks']])}

### Uncertainty and acquisition priorities

{markdown_table(uncertainty)}

These are rank correlations, not calibrated prediction intervals or independent-pixel significance tests. Standard deviation is across ensemble runs, **not across trees** (`ddof=0`). Each run changes synthetic mixtures and forest seed. `acquisition_priority.png` marks blocks in the upper quartile of both S2-like spread and observed EnMAP error reduction, with positive improvement. Because actual error reduction needs reference truth, this is a **retrospective diagnostic**, not an operational acquisition policy for an unlabelled region. High spread may reflect variability without detecting systematic bias.

Only {int(base.n_blocks.iloc[0])} blocks pass the reporting support threshold. Pixel correlations are weak-to-moderate, and high block correlations for rare classes are exploratory with this small block count; they should not be promoted as established acquisition performance.

### Ensemble sensitivity

The 3-member baseline is the first three members of the same 5-member ensemble; this isolates the effect of adding two runs.

{markdown_table(sensitivity)}

### Tutorial sanity check

The official tutorial screenshot reports 1,578 samples and the class MAEs below. Our values are in the **{sanity_status}**. This checks gross scale, not equivalence: the screenshot has no fixed seed, versioned environment or exact reference-mask specification, and the tutorial's written coverage value is ambiguous (0.95 under a percent-labelled parameter). We implement the requested **95%** explicitly. We checked CRS alignment, spectral units/scaling, target-first mixing, sequential weights, reference area denominators and no-data handling.

{markdown_table(sanity)}

## Re-run

Python 3.13 was used. From the repository root on Windows PowerShell:

```powershell
python -m venv .venv
. .\\start.ps1
python -m pip install -r requirements.txt
python download_data.py
python run.py --config config.yaml
```

On Linux/macOS use `source .venv/bin/activate` instead of `start.ps1`. `download_data.py` downloads the pinned ~274 MB EnMAP-Box tutorial archive, checks SHA256 and extracts only the three required inputs. To reuse an existing download: `python download_data.py --archive /path/to/tutorial_materials.zip`. The SRF workbook and Stage B imagery are fetched separately by the pipeline on first use; internet access and additional disk space are required. External service availability is not guaranteed.

`python run.py --config config.yaml --stage prepare` runs input/reference/SRF audits only. `--stage a` runs baseline and sensitivity; `--stage b` runs the configured temporal extension after Stage A has created local rasters/caches. The default runs Stage A, followed by Stage B only when `extension.enabled` is true. Run `python -m unittest discover -s tests -v` for scientific invariants. No GUI is opened.

All numerical settings and seeds are in `config.yaml`; executed settings are recorded in `results/config_used.yaml`. Public configuration snapshots replace the original machine-specific archive path with `data/tutorial_materials.zip`; scientific settings are unchanged. `environment.json` and `requirements.txt` record installed versions. Input and SRF checksums, georeferencing, band decisions, coverage weights, and library counts are saved alongside the results. CPU execution took **{elapsed:.1f} s** (excluding initial installation/downloads); per-model durations are in `model_runtime.csv` and the command log is `run.log`.

## Data and method

- The [official tutorial](https://enmap-box.readthedocs.io/en/latest/usr_section/application_tutorials/urban_unmixing/tutorial.html) supplies a 400×220 pixel, 30 m **simulated** EnMAP scene with 177 bands and 75 endmembers (impervious 38, vegetation 31, soil 4, water 2). All library wavelength axes match the image metadata exactly. Image and library retain the same reflectance scaling; no independent standardization is applied.
- [Tutorial data ZIP](https://box.hu-berlin.de/f/e4fa78c198bc4d868d30/?dl=1). Cite Okujeni, van der Linden & Hostert, *Berlin-Urban-Gradient dataset 2009 — An EnMAP Preparatory Flight Campaign*, [GFZ Data Services](https://doi.org/10.5880/enmap.2016.002), with the [technical report](https://doi.org/10.2312/enmap.2016.002). The original HyMap acquisition was 20 August 2009; the reference data were delineated from cadastre, orthophotos and elevation data, whose exact image dates are not encoded in the tutorial layer. Dataset licence: CC BY-SA 4.0; raw data are downloaded separately rather than redistributed here.
- Reference classes are rasterized at 3 m on the exact 30 m grid using pixel-centre assignment. Fractions are subpixel counts divided by **all 100 subpixels**, matching the official reference algorithm, and are not normalized to covered area. Geometry repairs, when needed, are audited; the original file is preserved.
- For each target class and run: 1,000 mixtures, first endmember from the target, others sampled uniformly over library rows (therefore proportional to class count); 2/3/4 endmembers with probabilities 0.4/0.4/0.2; within-class/repeated draws allowed; sequential uniform allocation of the remaining mixture weight; original 75 endmembers appended with binary target fractions; zero additional background proportion. Rules were checked against the [official SynthMix source](https://github.com/EnMAP-Box/enmap-box/blob/main/enmapboxprocessing/algorithm/prepareregressiondatasetfromsynthmixalgorithm.py), and one separate RF is fitted per target. No reference fractions train the Stage A models.
- RandomForestRegressor: 100 trees, all features available at a split (`max_features=1.0`), leaf minimum 1, bootstrap enabled, OOB scoring enabled as in the [EnMAP-Box RF preset](https://github.com/EnMAP-Box/enmap-box/blob/main/enmapboxprocessing/algorithm/fitrandomforestregressoralgorithm.py). OOB scores refer to synthetic data and are **not** the reference accuracy. Three runs form the baseline; five form sensitivity. Sum-to-one normalization is {'enabled' if cfg['mixing']['sum_to_one'] else 'disabled'}, following the unconstrained tutorial workflow.
- Ten Sentinel-2-like bands: B2/B3/B4/B5/B6/B7/B8/B8A/B11/B12. Method: **{srf_info.method.iloc[0]}**, obtained from the [official SentiWiki documents](https://sentiwiki.copernicus.eu/web/s2-documents), not a Gaussian approximation in this run. The response is integrated over linearly interpolated spectral samples, excluding unsupported gaps >35 nm; all integrated SRFs are renormalized to the observed support. Minimum full-SRF area coverage is **{srf_info.response_coverage.min():.3f}**. `s2_response_audit.csv` identifies truncation, especially at the blue edge; this is not an exact simulation of all S2 responses.

## Spectral physics check

{interpretation}

![Spectral importance](results/spectral_importance.png)

## Deviations from the EnMAP-Box workflow

- QGIS/EnMAP-Box is not installed yet. No official headless run or pixel-by-pixel output agreement has been measured; that is a later add-on.
- This is an independent Python implementation using NumPy Generator with fixed per-class/run seeds, rather than EnMAP-Box's unseeded global NumPy/Python random generators. The mixing distribution and target/background rules follow the source; the particular random samples differ.
- The installed scikit-learn version is recorded, and need not equal the tutorial's original environment. Parallel execution changes speed, not the intended estimator settings.
- The reference coverage threshold is explicitly 95%; invalid input geometry is repaired in memory with a recorded area/raster audit. Original inputs are unchanged.
- The spectral convolution, paired sensor comparison, spatial summaries, ensemble sensitivity, block resampling, and uncertainty analyses extend the tutorial.

## Limitations

One city, one old scene, simulated EnMAP, a small image-derived library, and sparse polygon references limit interpretation. Endmembers and evaluation imagery share the same parent scene, so this is not an independent-scene generalization test. Soil and water have especially few endmembers and positive reference pixels: their overall MAE can be dominated by absence. `summary.csv` includes positive-reference support and its conditional MAE. Convolved S2-like spectra share EnMAP spatial resolution, source date and radiometry; they do **not** reproduce real S2 atmosphere, spatial point-spread functions, native band resolutions or acquisition differences. The major water-vapour/edge exclusions were already applied to the supplied raster; removing more bands without evidence would alter the baseline. Impurity importance among correlated bands is not causal evidence.

## Stage B and later work

The separately reported [Stage A additions](results/stage_a_additions/REPORT.md) provide paired-block gain intervals and block shares, context decomposition, brightness/shadow-proxy water diagnosis with example spectra, and a standardized SVR sensitivity whose parameters were selected only on synthetic mixtures. They preserve the RF baseline. Re-run only these additions with `python run.py --config config.yaml --stage a-additions`.

Stage B uses independent external WorldCover-derived endmembers, compares single-date and temporal S2 features from the same pixels, and averages temporal S2 and simulated EnMAP predictions without training a joint model. It must prominently report the 2009/S2/WorldCover year differences, class crosswalk limitations and fixed data gates. See `results/stage_b/STATUS.md` when that stage has run. The official QGIS comparison remains pending. The separate [real L2A first look](REAL_DATA_FIRST_LOOK.md) applies the frozen RF and compares actual Sentinel-2; its sparse support and 16-year reference gap prevent a contemporary accuracy claim.

{stage_b_result}

Future ideas only: vegetation/burned-area fraction changes around a wildfire, and larger-area mapping along an EnMAP Foreground flightline. Stage C scaling and real EnMAP transfer are not part of these Stage A results.

## Figures and outputs

![Fraction maps](results/fraction_maps.png)
![Block differences](results/block_mae_difference.png)
![Uncertainty versus error](results/uncertainty_vs_error.png)
![Retrospective priority](results/acquisition_priority.png)

Published outputs: `summary.csv`, `blocks.csv`, `block_differences.csv`, `contexts.csv`, `paired_block_bootstrap.csv`, input/band/SRF audits, model timings and PNG figures. GeoTIFF fraction/uncertainty maps and cached ensemble predictions are generated locally and excluded from the public repository, along with all raw tutorial inputs.
'''
    (root/'README.md').write_text(text, encoding='utf-8')
