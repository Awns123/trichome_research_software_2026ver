from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

from common_utils import build_image_level_table

BASE = Path('/mnt/data')
OUT = BASE / 'figure2_trait_selection_validation_regenerated.png'

# Load image-level medians and reconstruct a size axis from size-dominant descriptors.
img = pd.read_csv(BASE / 'solanaceae_image_level_medians.csv')
img['curv_per_len'] = np.log1p(img['curvature_sum_rad'] / img['length_um'])

size_traits = [
    'area_um2', 'perimeter_um', 'ellipse_major_um', 'ellipse_minor_um',
    'length_um', 'straight_length_um', 'width_base_um', 'width_tip_um',
    'width_max_um', 'skeleton_length_um'
]
X_size = img[size_traits].astype(float).copy()
X_size = X_size.fillna(X_size.median(numeric_only=True))
X_size = (X_size - X_size.mean(axis=0)) / X_size.std(axis=0, ddof=0)
size_pc1 = PCA(n_components=1).fit_transform(X_size).ravel()

selected_traits = ['area_um2', 'skeleton_length_um', 'width_max_um', 'length_um', 'tortuosity', 'curvature_mean_rad', 'curv_per_len']
label_map = {
    'area_um2': 'Area',
    'skeleton_length_um': 'Skeleton Length',
    'width_max_um': 'Width (Max)',
    'length_um': 'Length',
    'tortuosity': 'Tortuosity',
    'curvature_mean_rad': 'Mean Curvature',
    'curv_per_len': 'Curvature per Length',
}
bar_rows = []
for t in selected_traits:
    rho, _ = spearmanr(img[t].astype(float), size_pc1)
    bar_rows.append({'trait': t, 'label': label_map[t], 'rho': float(rho)})
bar_df = pd.DataFrame(bar_rows)

# Correlation matrix for representative traits.
rep_traits = ['tortuosity', 'curvature_mean_rad', 'curv_per_len', 'area_um2', 'length_um', 'skeleton_length_um']
rep_corr = img[rep_traits].corr(method='spearman')
rep_corr.index = [label_map.get(x, x) for x in rep_corr.index]
rep_corr.columns = [label_map.get(x, x) for x in rep_corr.columns]

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

axes[0].barh(bar_df['label'], bar_df['rho'])
axes[0].axvline(0, linewidth=1)
axes[0].set_xlabel('Spearman correlation with size axis (PC1)')
axes[0].set_title('(A) Association with the size axis')

im = axes[1].imshow(rep_corr.values, vmin=-1, vmax=1)
axes[1].set_xticks(range(len(rep_corr.columns)))
axes[1].set_xticklabels(rep_corr.columns, rotation=45, ha='right')
axes[1].set_yticks(range(len(rep_corr.index)))
axes[1].set_yticklabels(rep_corr.index)
axes[1].set_title('(B) Redundancy among representative traits')
for i in range(rep_corr.shape[0]):
    for j in range(rep_corr.shape[1]):
        axes[1].text(j, i, f'{rep_corr.values[i, j]:.2f}', ha='center', va='center', fontsize=8)
fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)

fig.suptitle('Trait Selection Validation for the Focal Path–Curvature Block')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
