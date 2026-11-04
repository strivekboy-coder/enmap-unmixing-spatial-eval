# Stage B: fixed decisions before evaluating model errors

This is a method demonstration across a known time gap, not a same-date accuracy comparison: HyMap/simulated EnMAP is from 2009, candidate S2 seasons are 2017/2018, and WorldCover labels are from 2021. WorldCover agreement is a consistency screen, not evidence that a location remained unchanged throughout that interval.

## Fixed gates

1. External pure endmembers: >=30 per class, outside the entire tutorial raster footprint plus a 1 km exclusion, and within 10 km of that footprint. WorldCover is sampled onto a 10 m grid aligned to the tutorial 30 m grid; >=95% purity means all nine 10 m samples agree. No relaxed purity threshold or replacement using evaluation reference labels is allowed. At most 250 candidates/class are sampled with a fixed seed, and the >=30 gate applies again after cloud/feature filtering.
2. Reference consistency: WorldCover dominant mapped class agrees with the reference dominant class on >=70% of all 1,578 qualifying reference pixels; missing/unmapped WorldCover pixels count as disagreement. The stability-screen subset is exactly these agreeing pixels, defined before inspecting model error. Both full and screened results must be reported. This screen can preferentially select agreement with the external training labels and is therefore a sensitivity analysis, not a replacement for full-area results.
3. Clear observations: choose 2017 or 2018 by the larger median clear-observation count over the valid tutorial image footprint; break ties using the 10th percentile then earlier year. Require >=5 distinct clear dates on >=80% of valid image pixels. SCL classes 4/5/6 count as clear; snow, cloud, shadow and unclassified pixels do not. One source item per UTC calendar day is retained, preferring complete tutorial coverage, then lower scene cloud cover. Dates cannot be inflated by duplicate products/overlapping tiles.

The spatial region and numeric gates are in `config.yaml`. Gate failure is recorded and stops Stage B without manufacturing later results.

## Class crosswalk

WorldCover 50 built-up -> impervious; 60 bare/sparse vegetation -> soil; 80 permanent water -> water; 10/20/30/40/90/95/100 -> vegetation. Snow/ice and nodata are not endmembers. Built-up does not exactly match all impervious surfaces, and sparse vegetation is not necessarily pure bare soil. Map purity measures consistency of mapped labels, not physical subpixel purity.

## Models after gates pass

- (a) fixed Stage A EnMAP library model and (b) fixed Stage A S2-like convolved-library model: no tutorial reference training.
- (c0) best-date real S2 and (c) seasonal real S2: exactly the same external endmember pixels, target-first synthetic mixing, and **13-feature schema**. Best date is chosen by clear coverage among external candidate endmembers, independently of model error. Features are ten band medians and NDVI median/max/amplitude. For a one-date sequence, the band medians equal that date's bands, NDVI median=max and amplitude=0; for the full season they summarize all supported dates. This supplies the same engineered NDVI ratio to both models. Incomplete observations remain missing and cannot be imputed from evaluation reference fractions.
- (d) unweighted late fusion, `(a + c) / 2`, with no additional training. Ensemble-member fusion provides empirical spread; inverse-variance weighting is omitted unless separately justified, because a small spread can accompany systematic bias.
- Evaluate all models on matched supported tutorial pixels; report any extra cloud-driven exclusions and full versus fixed consistency-screen results. Save per-class, per-block and context MAE, ensemble-spread/error correlations, and block-wise `MAE(c) - MAE(d)` as a retrospective value-of-acquisition diagnostic.
- Monthly local clear fractions are estimated from SCL with a stated usable-area threshold. Expected attempts are `1/p` only under an independent-attempt, stationary-cloud approximation; report zero-success months without claiming a finite estimate. S2 cloudiness is a proxy for EnMAP conditions.

## Stop behaviour

Download or API errors are distinguished from scientific gate failures. Cached metadata allow retries. Stage A remains intact. Official QGIS comparison, Stage C scaling, and real EnMAP transfer remain later additions.

## Radiometry QC amendment, before any Stage B reference errors

Five selected historical images were reprocessed under baseline 05.00. Their `boa_offset_applied` flag is true, but `raster:bands` still advertises -0.1. Same-date comparisons against zero-offset legacy assets over 71,158 image pixels in red, NIR and SWIR show median raw-DN differences only -57 to +47; applying -0.1 introduces an artificial ~0.1 reflectance discontinuity. No reference labels or errors enter this check.

The five item IDs are explicitly listed in `config.yaml`. Effective offset zero is used only for these audited items, provided the provider flag remains true and the saved same-date audit includes all three bands with >=1,000 matched pixels and absolute median raw difference <=250 DN. All other items use their asset scale/offset. Both declared and effective offsets and source URLs are recorded. This fixes a demonstrated input-metadata inconsistency; it does not fit a radiometric correction to the evaluation targets. See `same_date_radiometric_audit.csv` and `check_radiometry.py`. The documented provider policy and a related metadata issue are [Earth Search notes](https://github.com/Element84/earth-search/blob/main/docs/collections/sentinel-2-l2a.md#gainoffset-in-items-after-jan-25-2022) and [issue 66](https://github.com/Element84/earth-search/issues/66).
