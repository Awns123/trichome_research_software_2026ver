from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common_utils import build_morph_distance_matrix, build_molecular_distance_matrix

BASE = Path('/mnt/data')
OUT = BASE / 'figure5_distance_heatmaps_regenerated.png'

morph = build_morph_distance_matrix(use_contraction=False)
mol = build_molecular_distance_matrix().loc[morph.index, morph.index]

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

im1 = axes[0].imshow(morph.values)
axes[0].set_xticks(range(len(morph.columns)))
axes[0].set_xticklabels(morph.columns, rotation=45, ha='right')
axes[0].set_yticks(range(len(morph.index)))
axes[0].set_yticklabels(morph.index)
axes[0].set_title('Morphological distance matrix')
fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)

im2 = axes[1].imshow(mol.values)
axes[1].set_xticks(range(len(mol.columns)))
axes[1].set_xticklabels(mol.columns, rotation=45, ha='right')
axes[1].set_yticks(range(len(mol.index)))
axes[1].set_yticklabels(mol.index)
axes[1].set_title('Molecular distance matrix')
fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)

fig.suptitle('Heatmap Representations of the Morphological and Molecular Distance Matrices')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
