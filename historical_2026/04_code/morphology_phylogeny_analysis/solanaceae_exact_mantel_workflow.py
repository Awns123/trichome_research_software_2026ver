from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr

BASE = Path('/mnt/data')
TREE_FILE = BASE / 'molecular_tree.txt'

RAW_SPECIES_FILES = {
    'B. suaveolens': BASE / 'B. suaveolens.csv',
    'C. annuum var. grossum': BASE / 'C. annuum var. grossum.csv',
    'L. chinense Mill.': BASE / 'L. chinense Mill..csv',
    "P. × hybrida ‘Dreams Red’": BASE / "P. × hybrida ‘Dreams Red’.csv",
    'P. alkekengi var. franchetii': BASE / 'P. alkekengi var. franchetii.csv',
    'S. nigrum': BASE / 'S. nigrum.csv',
}

TREE_LABEL_MAP = {
    'B. suaveolens': 'Brugmansiasuaveolens',
    'C. annuum var. grossum': 'CapsicumannuumL',
    'L. chinense Mill.': 'Lyciumchinense',
    "P. × hybrida ‘Dreams Red’": 'Petuniahybrida',
    'P. alkekengi var. franchetii': 'Physalisalkekengi',
    'S. nigrum': 'Solanumnigrum',
}


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    return (df - df.mean(axis=0)) / df.std(axis=0, ddof=0)


def upper_triangle_values(mat: np.ndarray) -> np.ndarray:
    idx = np.triu_indices_from(mat, k=1)
    return mat[idx]


def mantel_correlation(x_mat: np.ndarray, y_mat: np.ndarray, method: str = 'spearman') -> float:
    x = upper_triangle_values(x_mat)
    y = upper_triangle_values(y_mat)
    if method == 'spearman':
        return float(spearmanr(x, y).statistic)
    if method == 'pearson':
        return float(pearsonr(x, y).statistic)
    raise ValueError("method must be 'spearman' or 'pearson'")


def exact_mantel(x_mat: np.ndarray, y_mat: np.ndarray, method: str = 'spearman', alternative: str = 'greater'):
    obs = mantel_correlation(x_mat, y_mat, method=method)
    n = x_mat.shape[0]
    perm_stats = []
    for perm in itertools.permutations(range(n)):
        perm = np.array(perm)
        y_perm = y_mat[np.ix_(perm, perm)]
        perm_stats.append(mantel_correlation(x_mat, y_perm, method=method))
    perm_stats = np.asarray(perm_stats, dtype=float)
    if alternative == 'greater':
        p_value = float(np.mean(perm_stats >= obs - 1e-12))
    elif alternative == 'two-sided':
        p_value = float(np.mean(np.abs(perm_stats) >= abs(obs) - 1e-12))
    else:
        raise ValueError("alternative must be 'greater' or 'two-sided'")
    return obs, p_value, perm_stats


def build_image_level_table(use_contraction: bool = False) -> pd.DataFrame:
    rows = []
    for species_name, path in RAW_SPECIES_FILES.items():
        df = pd.read_csv(path)
        df['curv_per_len'] = df['curvature_sum_rad'] / df['length_um']
        img = df.groupby('image', dropna=False)[['tortuosity', 'curvature_mean_rad', 'curv_per_len', 'contraction_ratio']].median(numeric_only=True).copy()
        img['curv_per_len'] = np.log1p(img['curv_per_len'])
        if use_contraction:
            img['contraction_abs'] = np.log1p(np.abs(img['contraction_ratio']))
        img['species'] = species_name
        rows.append(img.reset_index())
    pooled = pd.concat(rows, ignore_index=True)
    trait_cols = ['tortuosity', 'curvature_mean_rad', 'curv_per_len']
    if use_contraction:
        trait_cols.append('contraction_abs')
    pooled[trait_cols] = zscore_df(pooled[trait_cols])
    return pooled


def build_species_centroids(image_level_table: pd.DataFrame, use_contraction: bool = False) -> pd.DataFrame:
    trait_cols = ['tortuosity', 'curvature_mean_rad', 'curv_per_len']
    if use_contraction:
        trait_cols.append('contraction_abs')
    return image_level_table.groupby('species')[trait_cols].median(numeric_only=True)


def build_distance_matrix(df: pd.DataFrame) -> pd.DataFrame:
    dist = squareform(pdist(df.values, metric='euclidean'))
    return pd.DataFrame(dist, index=df.index, columns=df.index)


def build_molecular_distance_matrix() -> pd.DataFrame:
    tree = Phylo.read(TREE_FILE, 'nexus')
    species_names = list(TREE_LABEL_MAP.keys())
    dist_mat = pd.DataFrame(index=species_names, columns=species_names, dtype=float)
    for sp_i in species_names:
        for sp_j in species_names:
            dist_mat.loc[sp_i, sp_j] = tree.distance(TREE_LABEL_MAP[sp_i], TREE_LABEL_MAP[sp_j])
    return dist_mat


if __name__ == '__main__':
    image_level = build_image_level_table(use_contraction=False)
    centroids = build_species_centroids(image_level, use_contraction=False)
    morph_dist = build_distance_matrix(centroids)
    mol_dist = build_molecular_distance_matrix().loc[centroids.index, centroids.index]

    s_r, s_p, _ = exact_mantel(morph_dist.values, mol_dist.values, method='spearman', alternative='greater')
    p_r, p_p, _ = exact_mantel(morph_dist.values, mol_dist.values, method='pearson', alternative='greater')
    print(f'Spearman exact Mantel: r_M = {s_r:.3f}, p = {s_p:.3f}')
    print(f'Pearson  exact Mantel: r_M = {p_r:.3f}, p = {p_p:.3f}')
