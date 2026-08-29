from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from common_utils import build_image_level_table, build_species_centroids, FOCAL_TRAITS

BASE = Path('/mnt/data')
OUT = BASE / 'figure3_pca_biplot_regenerated.png'

image_df = build_image_level_table(use_contraction=False)
centroids = build_species_centroids(image_df, use_contraction=False)

X = centroids[FOCAL_TRAITS].values
pca = PCA(n_components=2)
coords = pca.fit_transform(X)
loadings = pca.components_.T

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(coords[:, 0], coords[:, 1])
for i, species in enumerate(centroids.index):
    ax.text(coords[i, 0], coords[i, 1], species, fontsize=9)

scale = np.max(np.abs(coords)) * 0.8
for i, trait in enumerate(['Tortuosity', 'Mean Curvature', 'Curvature per Length']):
    ax.arrow(0, 0, loadings[i, 0] * scale, loadings[i, 1] * scale,
             head_width=0.05 * scale, length_includes_head=True)
    ax.text(loadings[i, 0] * scale * 1.08, loadings[i, 1] * scale * 1.08, trait, fontsize=9)

ax.axhline(0, linewidth=1)
ax.axvline(0, linewidth=1)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)')
ax.set_title('PCA Biplot of Species-Level Centroids')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
