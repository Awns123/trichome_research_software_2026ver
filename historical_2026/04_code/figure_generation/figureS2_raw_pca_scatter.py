from pathlib import Path

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from common_utils import build_image_level_table, FOCAL_TRAITS

BASE = Path('/mnt/data')
OUT = BASE / 'figureS2_raw_pca_scatter_regenerated.png'

image_df = build_image_level_table(use_contraction=False)
X = image_df[FOCAL_TRAITS].values
pca = PCA(n_components=2)
coords = pca.fit_transform(X)

fig, ax = plt.subplots(figsize=(8, 6))
for species, g in image_df.groupby('species'):
    idx = g.index.to_numpy()
    ax.scatter(coords[idx, 0], coords[idx, 1], label=species, s=18)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)')
ax.set_title('Raw Image-Level PCA Scatter Plot')
ax.legend(fontsize=7)
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
