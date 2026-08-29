from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path('/mnt/data')
OUT = BASE / 'figureS1_trait_correlation_heatmap_regenerated.png'

# Rebuild a compact but informative heatmap from the redundancy summary.
red_df = pd.read_csv(BASE / 'trait_redundancy_spearman_pairs.csv')
traits = sorted(set(red_df['trait1']).union(set(red_df['trait2'])))
mat = pd.DataFrame(np.eye(len(traits)), index=traits, columns=traits)
for _, row in red_df.iterrows():
    a, b, rho = row['trait1'], row['trait2'], row['spearman_rho']
    mat.loc[a, b] = rho
    mat.loc[b, a] = rho

fig, ax = plt.subplots(figsize=(12, 10))
im = ax.imshow(mat.values, vmin=-1, vmax=1)
ax.set_xticks(range(len(traits)))
ax.set_xticklabels(traits, rotation=90, fontsize=7)
ax.set_yticks(range(len(traits)))
ax.set_yticklabels(traits, fontsize=7)
ax.set_title('Supplementary Trait Correlation Heatmap')
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
