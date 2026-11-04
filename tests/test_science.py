"""Analytical invariants that catch scientifically consequential mistakes."""
import json
from pathlib import Path
import unittest

import numpy as np
import yaml

from data_io import response_weights
from evaluation import block_ids, rho
from mixing import mixtures

ROOT = Path(__file__).resolve().parents[1]


class ScientificInvariants(unittest.TestCase):
    def test_mixture_fraction_equals_one_hot_spectral_component(self):
        cfg = yaml.safe_load((ROOT/'config.yaml').read_text(encoding='utf-8'))['mixing']
        cfg['mixtures_per_class'] = 300
        # If each material has a unit spectrum, its abundance is observable
        # exactly in the corresponding spectral component, regardless of draws.
        library = np.eye(4, dtype='float32')
        labels = np.arange(4)
        for target in range(4):
            x, y = mixtures(library, labels, target, cfg, 42+target)
            np.testing.assert_allclose(x.sum(axis=1), 1, atol=1e-6)
            np.testing.assert_allclose(x[:, target], y, atol=1e-7)
            self.assertTrue(np.all(x >= 0))
            np.testing.assert_array_equal(x[-4:], library)
            np.testing.assert_array_equal(y[-4:], labels == target)

    def test_seed_reproducibility(self):
        cfg = yaml.safe_load((ROOT/'config.yaml').read_text(encoding='utf-8'))['mixing']
        cfg['mixtures_per_class'] = 20
        one = mixtures(np.eye(4), np.arange(4), 0, cfg, 44)
        two = mixtures(np.eye(4), np.arange(4), 0, cfg, 44)
        for a, b in zip(one, two):
            np.testing.assert_array_equal(a, b)

    def test_srf_preserves_constant_and_integrates_linear_spectrum(self):
        wavelengths = np.array([400., 410., 420., 430., 440.])
        axis = np.arange(405., 436.)
        response = np.ones((len(axis), 1))
        weights, coverage = response_weights(wavelengths, axis, response, 20)
        np.testing.assert_allclose(np.ones(5) @ weights, [1.], atol=1e-6)
        np.testing.assert_allclose(wavelengths @ weights, [420.], atol=1e-4)
        np.testing.assert_allclose(coverage, [1.])

    def test_srf_does_not_bridge_absorption_gap(self):
        wavelengths = np.array([400., 410., 490., 500.])
        axis = np.arange(400., 501.)
        responses = np.ones((len(axis), 1))
        weights, coverage = response_weights(wavelengths, axis, responses, 20)
        self.assertLess(coverage[0], .25)
        np.testing.assert_allclose(weights.sum(axis=0), [1.], atol=1e-6)

    def test_blocks_are_fixed_to_grid_and_edges_remain_partial(self):
        ids = block_ids((67, 67), 33)
        self.assertEqual(ids[0, 0], 0)
        self.assertEqual(ids[32, 32], 0)
        self.assertEqual(ids[33, 0], 3)
        self.assertEqual(ids[66, 66], 8)

    def test_spread_constant_is_not_perfect_confidence(self):
        self.assertTrue(np.isnan(rho(np.zeros(10), np.arange(10))))
        self.assertAlmostEqual(rho(np.arange(10), np.arange(10)), 1.)

    @unittest.skipUnless((ROOT/'results/cache/prepared.npz').exists(), 'Run prepare for input integration checks')
    def test_prepared_reference_has_consistent_masks_and_units(self):
        with np.load(ROOT/'results/cache/prepared.npz') as data:
            self.assertEqual(data['image'].shape, (177, 400, 220))
            self.assertEqual(data['library'].shape, (75, 177))
            self.assertTrue(np.isfinite(data['image'][:, data['valid_image']]).all())
            coverage = data['reference'].sum(axis=0)
            self.assertTrue(np.all(coverage <= 1.000001))
            self.assertTrue(np.all(coverage[data['valid']] >= .95-1e-7))
            self.assertTrue(np.all(data['valid_image'][data['valid']]))
            np.testing.assert_allclose(data['weights'].sum(axis=0), 1., atol=1e-6)
            # Image and spectral library have a common reflectance scale.
            ratio = np.nanmedian(data['image'][:, data['valid_image']]) / np.median(data['library'])
            self.assertTrue(.1 < ratio < 10)


if __name__ == '__main__':
    unittest.main()
