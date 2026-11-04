# Results notes

## Step A1–A2: inputs, reference grid and sensor response

- Original files remain unchanged; SHA-256 checksums are in `results/input_manifest.json`.
- 177/177 supplied bands retained; all 75 library wavelength axes match exactly.
- The two major water-vapour gaps were already removed from the tutorial input.
- 71,158 image pixels have finite values across all retained bands.
- Exactly **1,578** reference pixels meet 95% coverage, matching the official tutorial screenshot's sample count.
- Five self-intersecting geometries were repaired in memory. Relative area differences are < 3e-15, and **zero 3 m subpixels changed** in the rasterization comparison.
- ESA S2A SRF v4.0 was integrated for ten bands. All response-area coverages are at least **98.9708%**; the only meaningful truncation is B2 at the source blue edge.
- Analytical tests cover abundance conservation, deterministic mixing, linear/constant spectral response integration, missing-region handling, block alignment and prepared-data masks.

## Step A3–A7: baseline, sensitivity and spatial analysis

- Baseline ensemble=3; sensitivity=5 reuses the same first three members. Two sensors use paired synthetic mixtures and RF seeds. RF defaults follow the official 100-tree/OOB preset; no sum-to-one adjustment is applied.
- Baseline EnMAP MAE: impervious **0.11437**, vegetation **0.11129**, soil **0.09721**, water **0.24974**. Corresponding tutorial screenshot: 0.1179, 0.1088, 0.1094, 0.2514. Gross-scale reproduction is consistent; official program output comparison remains pending QGIS installation.
- Baseline S2-like macro MAE **0.15061**, EnMAP **0.14315**, improvement **0.00746** (0.746 percentage points). The largest macro gain by dominant class is in impervious-dominated locations (1.28 pp).
- All 12 blocks with any qualifying reference contribute to the paired bootstrap; 9 blocks have >=20 evaluation pixels and enter block maps/correlations.
- EnMAP macro pixel/block uncertainty–error Spearman rho: **0.315 / 0.583**. Only 9 supported blocks: correlations are exploratory.
- Water predictions remain badly biased high (~0.245 mean error), and EnMAP performs slightly worse than S2-like on water-positive pixels. Overall macro improvement must not hide this limitation.
- Stage A wall runtime: **146.3 seconds**, excluding installation/downloads. All three main figures plus spectral-importance and retrospective-priority figures were visually inspected.
- Stage A outputs and deviations are described in `README.md`; `results/summary.csv` is the machine-readable entry point.

## Stage A additions: gains, water diagnosis and SVR sensitivity

- The 2,000-resample paired-block macro gain is **0.007458**, 95% percentile CI **[0.004123, 0.009603]**. EnMAP is better in **11/12** blocks with any reference, or **8/9** supported blocks. Impervious-cover gain is unresolved by its interval **[-0.001452, 0.003648]**, with only 5/12 blocks showing improvement.
- Context decomposition includes dominant-reference class, impervious-gradient bins and vegetation-mixture bins. Overall soil improvement does not imply improvement on soil-dominated pixels, and overall water improvement does not imply improvement on water-positive pixels.
- Water diagnosis: overprediction increases with decreasing brightness, but 683 high false-water pixels lie outside the fixed dark proxy and 89 inside it. Two library water spectra and three deterministic example pixels are plotted. Target-first mixing supplies only 73 explicit zero-water training examples among 1,075 rows; its positive-water prevalence is reported as a plausible contributor, not a causal test. No correction or reference-based tuning was applied.
- Standardized RBF SVR used one common **C=100, epsilon=0.01, gamma=scale**, selected solely on independent synthetic mixtures. RF and SVR use the same three paired seeds; bounds [0,1] were predeclared for the SVR variant.
- SVR macro MAE: **EnMAP 0.108423**, **S2-like 0.126596**, gain **0.018174**. The average sign agrees with RF, while **impervious-cover gain flips sign** (EnMAP 0.107553 vs S2-like 0.104403). No universal class-wise advantage is claimed.
- All RF baseline output checksums are unchanged. Eight saved SVR class metrics, synthetic parameter ranking and both block-share denominators were independently checked; diagnosis figures were visually inspected. The supplementary report is `results/stage_a_additions/REPORT.md`.

## Stage B completed: independent temporal library and fixed late fusion

- 2018 was selected over 2017 by actual SCL observations: median clear dates 32 vs 8; 97.42% vs 70.06% of valid image pixels have >=5 clear dates. A season has 85 vs 25 distinct daily items after deduplication. The image-year gap is **2009 -> 2018 (9 years)**, and WorldCover labels are from **2021**.
- Independent endmembers outside the whole evaluation footprint plus 1 km, within 10 km: **239 impervious, 250 vegetation, 194 soil, 248 water** after common cloud/feature filtering. c0/c use the same pixels and 13-feature schema, with one-date NDVI median=max and amplitude=0. Best date: **2018-05-29**, selected from external-library cloud coverage only.
- WorldCover/reference dominant-class agreement is **81.43%** overall. **No soil-dominated reference pixel agrees with WorldCover** (0/39); this important class-wise limitation is explicit in the report. Consistency screening is not proof of stability.
- All five models use **1,495 common reference pixels** (83 exclusions); the consistency subset has **1,205**. Primary macro MAE: a EnMAP **0.1450**, b S2-like **0.1524**, c0 best-date S2 **0.1597**, c S2 temporal **0.1545**, d fixed 50:50 late fusion **0.1393**.
- Time-series gain is not uniform: impervious and vegetation worsen slightly relative to c0. Fusion improves macro error, but water worsens **0.0523 -> 0.1514** relative to c, and soil worsens **0.0995 -> 0.1773** relative to a. No reference-fitted weights or post-hoc correction replace these outcomes.
- Five reprocessed assets advertised an extra -0.1 offset despite already corrected DNs. Same-date red/NIR/SWIR comparisons without reference labels establish median raw differences -57 to +47 DN. Explicitly audited overrides prevent a false reflectance discontinuity; declared/effective offsets and source URLs remain recorded.
- Local strict SCL usability is 6/85 dates, giving a **14.17-attempt plug-in proxy** under strong independence/stationarity assumptions. Surface-dependent SCL behaviour and month sample counts limit interpretation; scene-level cloud-cover proxies are reported separately.
- The original download/processing attempt took **2,161 s**; cached processing took **~53–55 s**. Attempt history is retained rather than hiding the download cost. A subsequent cached Stage B command completed inside the restricted workspace with no additional network access.
- Twelve scientific invariant tests passed. Saved metrics, support masks, external geographic separation, radiometric rules and member-wise fusion were independently checked. Stage B figures were visually inspected. QGIS official output, real EnMAP and Stage C scaling remain later additions.

## Public repository preparation

- Added an MIT license for code and explicit data/result attribution. README now opens with the research question, measured findings and limitations; the report generator preserves this structure.
- Added `download_data.py`: pinned archive and per-input SHA256 verification, extraction of only three required tutorial files, preservation of existing verified inputs, and rejection of changed/corrupt files. The main configuration now uses a portable local archive path.
- Twelve scientific tests passed again; local-archive extraction, exact file selection, idempotency and corrupt-file rejection passed. Scientific models and reported measurements were not changed.
- Public release uses a fresh Git history containing code, documentation, configuration, tables, figures and small provenance records. Full local research history and rasters are preserved separately. Historical configuration snapshots in the public copy substitute only the machine-specific archive path; scientific settings remain unchanged.
