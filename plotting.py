"""Export figures with explicit masks, fraction units and shared colour scales."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from evaluation import block_ids, rho

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 170, 'figure.facecolor': 'white'})


def save(fig, path):
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def block_grid(table, field, shape, size):
    """Expand blocks onto the original grid so partial edge blocks stay honest."""
    data = np.full(shape, np.nan)
    ncols = int(np.ceil(shape[1]/size))
    for row in table.itertuples():
        r, c = divmod(row.block_id, ncols)
        data[r*size:min((r+1)*size, shape[0]), c*size:min((c+1)*size, shape[1])] = getattr(row, field)*100
    return data


def figures(cfg, out, reference, valid, profile, means, stds, wavelengths, importance):
    classes = cfg['classes']
    size = cfg['mixing']['ensemble_size']
    h, w = valid.shape
    extent = [0, w*.03, 0, h*.03]
    fig, axes = plt.subplots(3, 4, figsize=(12, 16), layout='constrained', sharex=True, sharey=True)
    for row, (name, data) in enumerate([('Reference (>=95% coverage)', np.where(valid[None], reference, np.nan)),
                                       ('EnMAP (simulated)', means['EnMAP']), ('S2-like (SRF convolution)', means['S2-like'])]):
        for col, class_name in enumerate(classes):
            ax = axes[row, col]
            im = ax.imshow(data[col], vmin=0, vmax=1, cmap='viridis', extent=extent, interpolation='nearest')
            ax.set_title(f'{class_name.capitalize()}' if row == 0 else '')
            ax.set_facecolor('#eeeeee')
            if col == 0:
                ax.set_ylabel(name+'\nNorthing from southern grid edge (km)')
            if row == 2:
                ax.set_xlabel('Easting from grid origin (km)')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=.5, label='Land-cover fraction (0-1)')
    fig.suptitle(f'Berlin tutorial | 30 m grid | ensemble {size}\nGrey = no image data / no qualifying reference', fontsize=15)
    save(fig, out/'fraction_maps.png')

    blocks = pd.read_csv(out/'block_differences.csv')
    blocks = blocks[blocks.ensemble_size == size]
    block_size = cfg['evaluation']['block_size_pixels']
    rows, cols = int(np.ceil(h/block_size)), int(np.ceil(w/block_size))
    max_abs = max(.01, float(blocks.s2_minus_enmap_mae.abs().max()*100))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)
    fig, axes = plt.subplots(1, 5, figsize=(14, 6), layout='constrained', sharex=True, sharey=True)
    for ax, name in zip(axes, classes+['macro']):
        selected = blocks[blocks['class'] == name]
        data = block_grid(selected, 's2_minus_enmap_mae', (h, w), block_size)
        im = ax.imshow(data, cmap='RdBu', norm=norm, extent=extent, interpolation='nearest')
        ax.set_title(name.capitalize())
        ax.set_facecolor('#eeeeee')
        ax.set_xlabel('Easting (km)')
    axes[0].set_ylabel('Northing (km)')
    fig.colorbar(im, ax=axes.tolist(), shrink=.8, label='S2-like MAE minus EnMAP MAE (percentage points)')
    fig.suptitle(f'Where hyperspectral information helps | 33 x 33 pixel blocks\nBlue: lower EnMAP error; red: lower S2-like error; at least {cfg["evaluation"]["minimum_pixels_per_block"]} reference pixels/block', fontsize=13)
    save(fig, out/'block_mae_difference.png')

    fig, axes = plt.subplots(3, 4, figsize=(14, 10), layout='constrained')
    block_table = pd.read_csv(out/'blocks.csv')
    block_table = block_table[block_table.ensemble_size == size]
    rng = np.random.default_rng(cfg['seed'])
    indices = rng.choice(valid.sum(), size=min(valid.sum(), cfg['evaluation']['plot_sample_pixels']), replace=False)
    for col, name in enumerate(classes):
        for row, sensor in enumerate(means):
            unc = stds[sensor][col, valid]
            err = np.abs(means[sensor][col, valid]-reference[col, valid])
            ax = axes[row, col]
            ax.hexbin(unc[indices]*100, err[indices]*100, gridsize=26, mincnt=1, bins='log', cmap='Blues')
            ax.set_title(f'{sensor} | {name}\npixel Spearman rho={rho(unc, err):.2f}')
            ax.set_xlabel('Ensemble std (pp)'); ax.set_ylabel('Absolute error (pp)')
        ax = axes[2, col]
        for sensor, colour in [('EnMAP', '#146c94'), ('S2-like', '#d56b25')]:
            table = block_table[(block_table.sensor == sensor) & (block_table['class'] == name)]
            ax.scatter(table.mean_uncertainty*100, table.mae*100, color=colour,
                       label=f'{sensor}: rho={rho(table.mean_uncertainty, table.mae):.2f}', s=28, alpha=.85)
        ax.set_title(f'{name} | block means'); ax.set_xlabel('Block mean std (pp)'); ax.set_ylabel('Block MAE (pp)')
        ax.legend(fontsize=8)
    fig.suptitle('Does ensemble spread locate prediction errors?\nSpread reflects synthetic-training/model variability, not calibrated confidence', fontsize=14)
    save(fig, out/'uncertainty_vs_error.png')

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), layout='constrained', sharex=True)
    breaks = np.flatnonzero(np.diff(wavelengths) > cfg['spectral']['maximum_interpolation_gap_nm']) + 1
    for ax, target in zip(axes, [0, 2]):
        for region, colour, label in [( (1350, 1450), '#aaa', 'Water vapour'), ((1800, 1950), '#aaa', None),
                                      ((690, 750), '#a9cbae', 'Red edge'), ((2180, 2220), '#e9c982', '~2200 nm'),
                                      ((2280, 2320), '#dca9a1', '~2300 nm')]:
            ax.axvspan(*region, color=colour, alpha=.3, label=label)
        for section in np.split(np.arange(len(wavelengths)), breaks):
            ax.plot(wavelengths[section], importance[target, section], color='#205b78', linewidth=1.4)
        ax.set_ylabel('RF impurity importance')
        ax.set_title(f'{classes[target].capitalize()} fraction | correlated wavelengths limit physical attribution')
        ax.legend(loc='upper left', ncol=4, fontsize=8)
    axes[-1].set_xlabel('Wavelength (nm)')
    fig.suptitle('Spectral physics audit | gaps were already absent from the tutorial raster', fontsize=14)
    save(fig, out/'spectral_importance.png')

    selected = blocks[blocks['class'] == 'macro'].copy()
    fig, axes = plt.subplots(1, 2, figsize=(8, 7), layout='constrained')
    for ax, field, title, cmap in [(axes[0], 'mean_uncertainty_s2', 'S2-like ensemble spread', 'magma'),
                                    (axes[1], 's2_minus_enmap_mae', 'Observed EnMAP error reduction', 'RdBu')]:
        data = block_grid(selected, field, (h, w), block_size)
        options = {'norm': norm} if field == 's2_minus_enmap_mae' else {}
        im = ax.imshow(data, extent=extent, cmap=cmap, interpolation='nearest', **options)
        ax.set_title(title); ax.set_xlabel('Easting (km)'); ax.set_facecolor('#eeeeee')
        for row in selected[selected.priority_retrospective].itertuples():
            r, c = divmod(row.block_id, cols)
            center_x = (c*block_size + min((c+1)*block_size, w))*.03/2
            center_y = extent[3]-(r*block_size+min((r+1)*block_size, h))*.03/2
            ax.plot(center_x, center_y, '*', color='#36bd52', markersize=13, markeredgecolor='black')
        fig.colorbar(im, ax=ax, shrink=.7, label='Percentage points')
    axes[0].set_ylabel('Northing (km)')
    fig.suptitle('Retrospective acquisition-priority diagnostic\nStars: upper-quartile spread AND error reduction; requires reference truth', fontsize=12)
    save(fig, out/'acquisition_priority.png')
