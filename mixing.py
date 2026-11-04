"""Target-first synthetic mixtures following EnMAP-Box's sequential weights."""
from __future__ import annotations

import logging
import time
import numpy as np
from sklearn.ensemble import RandomForestRegressor

LOG = logging.getLogger(__name__)


def mixtures(library, labels, target, settings, seed):
    if settings['background_proportion'] != 0:
        raise ValueError('Only the tutorial default (zero added background mixtures) is supported.')
    rng = np.random.default_rng(seed)
    count = settings['mixtures_per_class']
    xs, ys = [], []
    target_indices = np.flatnonzero(labels == target)
    if not len(target_indices):
        raise ValueError(f'No endmembers for target class {target}.')
    # Uniform sampling of all library rows is exactly proportional to class size.
    for _ in range(count):
        complexity = rng.choice(settings['complexities'], p=settings['probabilities'])
        indices = np.concatenate(([rng.choice(target_indices)], rng.integers(0, len(labels), complexity-1)))
        weights = []
        remainder = 1.
        for _ in range(complexity-1):
            weight = rng.random() * remainder
            weights.append(weight)
            remainder -= weight
        weights.append(remainder)
        weights = np.asarray(weights)
        xs.append(weights @ library[indices])
        ys.append(weights @ (labels[indices] == target))
    if settings['include_original_endmembers']:
        xs.extend(library)
        ys.extend((labels == target).astype(float))
    return np.asarray(xs, dtype='float32'), np.asarray(ys, dtype='float32')


def train_sensor_pair(cfg, image, library, labels, valid_image, weights, out):
    """Reuse exactly the same mixtures and RF seeds for both sensor variants."""
    baseline = cfg['mixing']['ensemble_size']
    total = cfg['mixing']['sensitivity_ensemble_size']
    nclasses = len(cfg['classes'])
    nvalid = int(valid_image.sum())
    pixels = image[:, valid_image].T
    sensor_pixels = {'EnMAP': pixels, 'S2-like': pixels @ weights}
    predictions = {s: np.zeros((total, nclasses, nvalid), dtype='float32') for s in sensor_pixels}
    importance = {s: np.zeros((total, nclasses, x.shape[1])) for s, x in sensor_pixels.items()}
    timings = []
    for run in range(total):
        for target, name in enumerate(cfg['classes']):
            seed = cfg['seed'] + run*100 + target
            x, y = mixtures(library, labels, target, cfg['mixing'], seed)
            for sensor in sensor_pixels:
                tic = time.perf_counter()
                train_x = x if sensor == 'EnMAP' else x @ weights
                model = RandomForestRegressor(**cfg['forest'], random_state=seed)
                model.fit(train_x, y)
                predictions[sensor][run, target] = model.predict(sensor_pixels[sensor])
                importance[sensor][run, target] = model.feature_importances_
                elapsed = time.perf_counter()-tic
                timings.append({'sensor': sensor, 'run': run+1, 'class': name,
                                'seconds': elapsed, 'oob_r2_synthetic': model.oob_score_, 'seed': seed})
                LOG.info('run %d/%d %s %s: %.2fs, synthetic OOB R2=%.3f', run+1, total, sensor, name, elapsed, model.oob_score_)
        np.savez_compressed(out / 'cache' / 'ensemble_predictions.npz', **predictions)
    # If enabled, normalize each member so that its mean and std are coherent.
    if cfg['mixing']['sum_to_one']:
        for sensor in predictions:
            denom = predictions[sensor].sum(axis=1, keepdims=True)
            predictions[sensor] /= np.where(denom > 0, denom, 1.)
    return predictions, importance, timings
