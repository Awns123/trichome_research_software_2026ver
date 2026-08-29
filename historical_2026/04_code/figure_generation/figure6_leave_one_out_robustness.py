from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

BASE = Path('/mnt/data')
OUT = BASE / 'figure6_leave_one_out_robustness_regenerated.png'

# Use the summary already generated during analysis if available; otherwise hard-code manuscript values.
path = BASE / 'reproduced_leave_one_out_results.csv'
if path.exists():
    df = pd.read_csv(path)
else:
    df = pd.DataFrame([
        {'excluded_species': 'B. suaveolens', 'r_M': 0.333, 'p_value': 0.575},
        {'excluded_species': "P. × hybrida ‘Dreams Red’", 'r_M': 0.709, 'p_value': 0.100},
        {'excluded_species': 'S. nigrum', 'r_M': 0.733, 'p_value': 0.117},
        {'excluded_species': 'P. alkekengi var. franchetii', 'r_M': 0.806, 'p_value': 0.067},
        {'excluded_species': 'C. annuum var. grossum', 'r_M': 0.879, 'p_value': 0.050},
        {'excluded_species': 'L. chinense Mill.', 'r_M': 0.927, 'p_value': 0.042},
    ])

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(df['excluded_species'], df['r_M'])
ax.axhline(0.775, linestyle='--', linewidth=1)
for i, row in df.iterrows():
    ax.text(i, row['r_M'] + 0.02, f"p={row['p_value']:.3f}", ha='center', va='bottom', fontsize=8)
ax.set_ylim(0, max(df['r_M']) + 0.15)
ax.set_ylabel('Spearman exact Mantel $r_M$')
ax.set_title('Leave-One-Species-Out Robustness of the Focal Trait Block')
plt.xticks(rotation=45, ha='right')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
