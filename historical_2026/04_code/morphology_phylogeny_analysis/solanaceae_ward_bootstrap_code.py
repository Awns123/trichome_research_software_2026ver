from __future__ import annotations

from collections import Counter
from pathlib import Path
import itertools

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist

BASE = Path('/mnt/data')

RAW_SPECIES_FILES = {
    'B. suaveolens': BASE / 'B. suaveolens.csv',
    'C. annuum var. grossum': BASE / 'C. annuum var. grossum.csv',
    'L. chinense Mill.': BASE / 'L. chinense Mill..csv',
    "P. × hybrida ‘Dreams Red’": BASE / "P. × hybrida ‘Dreams Red’.csv",
    'P. alkekengi var. franchetii': BASE / 'P. alkekengi var. franchetii.csv',
    'S. nigrum': BASE / 'S. nigrum.csv',
}


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    return (df - df.mean(axis=0)) / df.std(axis=0, ddof=0)


def build_image_level_table() -> pd.DataFrame:
    rows = []
    for species_name, path in RAW_SPECIES_FILES.items():
        df = pd.read_csv(path)
        df['curv_per_len'] = np.log1p(df['curvature_sum_rad'] / df['length_um'])
        img = df.groupby('image', dropna=False)[['tortuosity', 'curvature_mean_rad', 'curv_per_len']].median(numeric_only=True).copy()
        img['species'] = species_name
        rows.append(img.reset_index())
    pooled = pd.concat(rows, ignore_index=True)
    pooled[['tortuosity', 'curvature_mean_rad', 'curv_per_len']] = zscore_df(
        pooled[['tortuosity', 'curvature_mean_rad', 'curv_per_len']]
    )
    return pooled


def species_centroids_from_bootstrap(image_df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    parts = []
    for sp, g in image_df.groupby('species'):
        idx = rng.choice(g.index.to_numpy(), size=len(g), replace=True)
        g_boot = g.loc[idx, ['tortuosity', 'curvature_mean_rad', 'curv_per_len']]
        parts.append((sp, g_boot.median(axis=0)))
    cent = pd.DataFrame({sp: vals for sp, vals in parts}).T
    cent.index.name = 'species'
    return cent


def tree_clades_from_linkage(Z, labels):
    root, nodes = to_tree(Z, rd=True)
    id_to_labels = {i: frozenset([labels[i]]) for i in range(len(labels))}
    clades = []
    for node in nodes:
        if node.id < len(labels):
            continue
        left = id_to_labels[node.left.id]
        right = id_to_labels[node.right.id]
        merged = left | right
        id_to_labels[node.id] = merged
        if 1 < len(merged) < len(labels):
            clades.append(merged)
    return frozenset(clades)


if __name__ == '__main__':
    image_df = build_image_level_table()
    labels = list(RAW_SPECIES_FILES.keys())

    # original topology
    orig_cent = image_df.groupby('species')[['tortuosity', 'curvature_mean_rad', 'curv_per_len']].median(numeric_only=True).loc[labels]
    Z_orig = linkage(pdist(orig_cent.values, metric='euclidean'), method='ward')
    orig_clades = tree_clades_from_linkage(Z_orig, labels)

    rng = np.random.default_rng(123)
    n_boot = 10000
    counts = Counter()
    topologies = Counter()

    for _ in range(n_boot):
        cent = species_centroids_from_bootstrap(image_df, rng).loc[labels]
        Z = linkage(pdist(cent.values, metric='euclidean'), method='ward')
        clades = tree_clades_from_linkage(Z, labels)
        topologies[clades] += 1
        for cl in orig_clades:
            if cl in clades:
                counts[cl] += 1

    for cl in sorted(orig_clades, key=lambda x: (len(x), sorted(x))):
        print(f"{' | '.join(sorted(cl))}: {counts[cl] / n_boot * 100:.2f}%")
