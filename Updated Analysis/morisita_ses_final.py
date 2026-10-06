"""
morisita_ses_097.py

Standardised effect size (SES) of the pairwise Morisita index for every
cell-type pair, using the two Morisita implementations in this project.

Output 1: <prefix>_SES_grid_shuffle.csv   (analyze_cell_pairs_grid_shuffle...py)
    - continuous Morisita-Horn formula, binning min(int(frac*g), g-1)
    - null: each cell type's grid-count column shuffled independently,
      using EXACTLY the same RNG calls as pairwise_morisita_index_test()
      (default_rng(42), rng.shuffle per column per permutation), so
      morisita_horn / p_value / fdr / significant reproduce that script.

Output 2: <prefix>_SES_label_shuffle.csv  (demo_setup_leiden_sp_data.py)
    - Simpson-corrected Morisita, binning int(frac*(g-1))
    - null: cell-label shuffle (coordinates fixed, leiden labels permuted)

Columns (both files):
    cat1, cat2, morisita_horn | morisita, p_value, fdr, significant,
    null_mean, null_sd, ses, p_two_sided, fdr_two_sided

    p_value       one-sided (null >= observed), (count + 1) / (n_perm + 1)
    fdr           Benjamini-Hochberg on p_value
    significant   fdr < ALPHA
    ses           (observed - null_mean) / null_sd   (null_sd with ddof=1)
    p_two_sided   min(1, 2 * min(p_greater, p_less))
    fdr_two_sided Benjamini-Hochberg on p_two_sided

Usage:
    python morisita_ses_097.py 097metadataUniform32125.csv
    python morisita_ses_097.py <file.csv> <output_prefix>
"""

import os
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

GRID_SIZE = 10
N_PERM = 1000
SEED = 42
ALPHA = 0.05


# ---------------------------------------------------------------- binning
def box_ids_label_method(x, y, g):
    """Binning of calculate_morisita_with_sortcategory() (demo_setup)."""
    xb = ((x - x.min()) / (x.max() - x.min()) * (g - 1)).astype(int)
    yb = ((y - y.min()) / (y.max() - y.min()) * (g - 1)).astype(int)
    return xb * g + yb


def box_ids_grid_method(x, y, g):
    """Binning of build_abundance_matrix() (grid_shuffle script)."""
    xb = np.minimum(((x - x.min()) / (x.max() - x.min()) * g).astype(int), g - 1)
    yb = np.minimum(((y - y.min()) / (y.max() - y.min()) * g).astype(int), g - 1)
    return xb * g + yb


def count_matrix(box_id, labels, n_boxes, k):
    """(n_boxes, k) matrix of cell counts per grid box per cell type."""
    return np.bincount(box_id * k + labels, minlength=n_boxes * k).reshape(n_boxes, k).astype(float)


# ---------------------------------------------------------------- formulas
def morisita_simpson_all(C):
    """Simpson-corrected Morisita for all pairs. C: (boxes, k) or (perm, boxes, k)."""
    n = C.sum(axis=-2)
    simpson = (C * (C - 1)).sum(axis=-2) / (n * (n - 1))
    cross = np.einsum('...bi,...bj->...ij', C, C)
    return 2 * cross / ((simpson[..., :, None] + simpson[..., None, :]) * n[..., :, None] * n[..., None, :])


def morisita_horn_all(C):
    """Continuous Morisita-Horn for all pairs (formula of morisita_index())."""
    n = C.sum(axis=-2)
    d = (C ** 2).sum(axis=-2) / n ** 2
    cross = np.einsum('...bi,...bj->...ij', C, C)
    return 2 * cross / ((d[..., :, None] + d[..., None, :]) * n[..., :, None] * n[..., None, :])


# ---------------------------------------------------------------- nulls
def null_grid_shuffle(C, n_perm, seed):
    """Same shuffle sequence as pairwise_morisita_index_test()."""
    rng = np.random.default_rng(seed)
    perms = np.empty((n_perm,) + C.shape)
    for p in range(n_perm):
        P = C.copy()
        for col in range(C.shape[1]):
            rng.shuffle(P[:, col])
        perms[p] = P
    return morisita_horn_all(perms)                            # (n_perm, k, k)


def null_label_shuffle(box_id, labels, n_boxes, k, n_perm, seed, batch=50):
    rng = np.random.default_rng(seed)
    out = []
    for start in range(0, n_perm, batch):
        mats = [count_matrix(box_id, rng.permutation(labels), n_boxes, k)
                for _ in range(min(batch, n_perm - start))]
        out.append(morisita_simpson_all(np.stack(mats)))
    return np.concatenate(out)


# ---------------------------------------------------------------- summary
def summarise(obs, null, cats, value_col):
    rows = []
    n_perm = null.shape[0]
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            nv, o = null[:, i, j], obs[i, j]
            mu, sd = nv.mean(), nv.std(ddof=1)
            p_greater = (np.sum(nv >= o) + 1) / (n_perm + 1)
            p_less = (np.sum(nv <= o) + 1) / (n_perm + 1)
            rows.append({
                'cat1': cats[i], 'cat2': cats[j], value_col: o,
                'p_value': p_greater,
                'null_mean': mu, 'null_sd': sd,
                'ses': (o - mu) / sd if sd > 0 else np.nan,
                'p_two_sided': min(1.0, 2 * min(p_greater, p_less)),
            })
    df = pd.DataFrame(rows)
    reject, fdr, _, _ = multipletests(df['p_value'], method='fdr_bh', alpha=ALPHA)
    df['fdr'] = fdr
    df['significant'] = reject
    df['fdr_two_sided'] = multipletests(df['p_two_sided'], method='fdr_bh')[1]
    return df[['cat1', 'cat2', value_col, 'p_value', 'fdr', 'significant',
               'null_mean', 'null_sd', 'ses', 'p_two_sided', 'fdr_two_sided']]


def main(path, prefix):
    df = pd.read_csv(path)
    cats = sorted(df['leiden'].unique().tolist())              # same order as the original script
    k = len(cats)
    lab = df['leiden'].map({c: i for i, c in enumerate(cats)}).values
    x, y = df['x_centroid'].values, df['y_centroid'].values
    nb = GRID_SIZE * GRID_SIZE

    # grid-shuffle script (Morisita-Horn)
    C = count_matrix(box_ids_grid_method(x, y, GRID_SIZE), lab, nb, k)
    grid = summarise(morisita_horn_all(C), null_grid_shuffle(C, N_PERM, SEED), cats, 'morisita_horn')
    grid.to_csv(f'{prefix}_SES_grid_shuffle.csv', index=False)

    # demo_setup (Simpson-corrected Morisita)
    bid = box_ids_label_method(x, y, GRID_SIZE)
    obs = morisita_simpson_all(count_matrix(bid, lab, nb, k))
    label = summarise(obs, null_label_shuffle(bid, lab, nb, k, N_PERM, SEED), cats, 'morisita')
    label.to_csv(f'{prefix}_SES_label_shuffle.csv', index=False)

    return grid, label


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'Metadata_S4_TOP_Uniform_23025.csv'
    prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(os.path.basename(path))[0]
    grid, label = main(path, prefix)
    pd.set_option('display.width', 250)
    print(f'--- {prefix}_SES_grid_shuffle.csv ---')
    print(grid.to_string(index=False))
    print(f'\n--- {prefix}_SES_label_shuffle.csv ---')
    print(label.to_string(index=False))
