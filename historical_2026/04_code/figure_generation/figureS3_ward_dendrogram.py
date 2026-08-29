from pathlib import Path

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from common_utils import build_image_level_table, build_species_centroids, FOCAL_TRAITS

BASE = Path('/mnt/data')
OUT = BASE / 'figureS3_ward_dendrogram_regenerated.png'

image_df = build_image_level_table(use_contraction=False)
centroids = build_species_centroids(image_df, use_contraction=False)
Z = linkage(pdist(centroids[FOCAL_TRAITS].values, metric='euclidean'), method='ward')

fig, ax = plt.subplots(figsize=(10, 5))
dendrogram(Z, labels=centroids.index.tolist(), ax=ax, orientation='top')
ax.set_ylabel('Distance')
ax.set_title('Ward Hierarchical Clustering Dendrogram')
plt.xticks(rotation=45, ha='right')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
