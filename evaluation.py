"""Matched-pixel metrics, spatial contexts, and paired block resampling."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score


def rho(x, y):
    x, y = np.asarray(x), np.asarray(y)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def block_ids(shape, size):
    rows, cols = np.indices(shape)
    return (rows // size) * int(np.ceil(shape[1]/size)) + cols // size


def summarize(cfg, reference, valid, means, stds, ensemble_size):
    classes = cfg['classes']
    ids = block_ids(valid.shape, cfg['evaluation']['block_size_pixels'])[valid]
    ref = reference[:, valid]
    dominant = ref.argmax(axis=0)
    edges = np.asarray(cfg['evaluation']['gradient_edges'])
    gradient = np.clip(np.digitize(ref[0], edges)-1, 0, len(edges)-2)
    summary, blocks, contexts = [], [], []
    for sensor in means:
        pred, uncertainty = means[sensor][:, valid], stds[sensor][:, valid]
        abs_error = np.abs(pred-ref)
        for target, name in list(enumerate(classes)) + [(None, 'macro')]:
            err = abs_error.mean(axis=0) if target is None else abs_error[target]
            unc = uncertainty.mean(axis=0) if target is None else uncertainty[target]
            overall = float(err.mean())
            block_rows = []
            for bid in np.unique(ids):
                mask = ids == bid
                if mask.sum() < cfg['evaluation']['minimum_pixels_per_block']:
                    continue
                b = {'sensor': sensor, 'ensemble_size': ensemble_size, 'class': name,
                     'block_id': int(bid), 'n_pixels': int(mask.sum()), 'mae': float(err[mask].mean()),
                     'mean_uncertainty': float(unc[mask].mean()),
                     'mean_reference_impervious': float(ref[0, mask].mean())}
                blocks.append(b)
                block_rows.append(b)
            bm = np.array([b['mae'] for b in block_rows])
            bu = np.array([b['mean_uncertainty'] for b in block_rows])
            row = {'sensor': sensor, 'ensemble_size': ensemble_size, 'class': name,
                   'n_pixels': int(valid.sum()), 'mae': overall,
                   'rmse': float(np.sqrt(np.mean((pred-ref)**2))) if target is None else float(np.sqrt(np.mean((pred[target]-ref[target])**2))),
                   'bias': float((pred-ref).mean()) if target is None else float((pred[target]-ref[target]).mean()),
                   'r2': np.nan if target is None else r2_score(ref[target], pred[target]),
                   'reference_positive_pixels': np.nan if target is None else int(np.sum(ref[target] > 0)),
                   'mae_reference_positive': np.nan if target is None or not np.any(ref[target] > 0) else float(err[ref[target] > 0].mean()),
                   'n_blocks': len(bm), 'block_mae_min': bm.min() if len(bm) else np.nan,
                   'block_mae_median': np.median(bm) if len(bm) else np.nan,
                   'block_mae_max': bm.max() if len(bm) else np.nan,
                   'share_blocks_above_twice_overall': float((bm > 2*overall).mean()) if len(bm) else np.nan,
                   'spearman_pixel_std_abs_error': rho(unc, err),
                   'spearman_block_std_mae': rho(bu, bm)}
            summary.append(row)
            groupings = [('dominant_reference', dominant, classes),
                         ('reference_impervious_fraction', gradient, ['0-25%', '25-50%', '50-75%', '75-100%'])]
            for kind, groups, labels in groupings:
                for group, label in enumerate(labels):
                    mask = groups == group
                    if mask.any():
                        contexts.append({'sensor': sensor, 'ensemble_size': ensemble_size, 'class': name,
                                         'context_type': kind, 'context': label, 'n_pixels': int(mask.sum()),
                                         'mae': float(err[mask].mean())})
    return pd.DataFrame(summary), pd.DataFrame(blocks), pd.DataFrame(contexts)


def differences(cfg, reference, valid, means, blocks, ensemble_size):
    a = blocks[blocks.sensor == 'EnMAP'].drop(columns=['sensor']).copy()
    b = blocks[blocks.sensor == 'S2-like'].drop(columns=['sensor']).copy()
    keys = ['ensemble_size', 'class', 'block_id']
    diff = a.merge(b, on=keys, suffixes=('_enmap', '_s2'), validate='one_to_one')
    diff['s2_minus_enmap_mae'] = diff.mae_s2 - diff.mae_enmap
    q = cfg['evaluation']['priority_quantile']
    # Retrospective diagnostic: error reduction requires reference data and is
    # not a deployable acquisition policy in an unlabelled region.
    diff['priority_retrospective'] = False
    for name in diff['class'].unique():
        mask = diff['class'] == name
        rows = diff.loc[mask]
        selected = ((rows.mean_uncertainty_s2 >= rows.mean_uncertainty_s2.quantile(q)) &
                    (rows.s2_minus_enmap_mae >= rows.s2_minus_enmap_mae.quantile(q)) &
                    (rows.s2_minus_enmap_mae > 0))
        diff.loc[rows.index, 'priority_retrospective'] = selected
    # Paired block bootstrap uses all evaluation pixels, including small blocks;
    # plotting/reporting block thresholds do not redefine the overall metric.
    ids = block_ids(valid.shape, cfg['evaluation']['block_size_pixels'])[valid]
    ae = np.abs(means['EnMAP'][:, valid] - reference[:, valid])
    be = np.abs(means['S2-like'][:, valid] - reference[:, valid])
    rng = np.random.default_rng(cfg['seed'])
    unique, inverse = np.unique(ids, return_inverse=True)
    counts = np.bincount(inverse)
    draw = rng.integers(0, len(unique), size=(cfg['evaluation']['bootstrap_resamples'], len(unique)))
    intervals = []
    for target, name in list(enumerate(cfg['classes'])) + [(None, 'macro')]:
        delta = (be-ae).mean(axis=0) if target is None else be[target]-ae[target]
        sums = np.bincount(inverse, weights=delta)
        sample = sums[draw].sum(axis=1) / counts[draw].sum(axis=1)
        intervals.append({'ensemble_size': ensemble_size, 'class': name,
                          's2_minus_enmap_mae': float(delta.mean()), 'ci025': np.quantile(sample, .025),
                          'ci975': np.quantile(sample, .975), 'n_spatial_blocks': len(unique)})
    return diff, pd.DataFrame(intervals)
