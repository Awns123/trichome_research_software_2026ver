from pathlib import Path

import matplotlib.pyplot as plt
from Bio import Phylo

from common_utils import TREE_FILE, RAW_SPECIES_FILES, TREE_LABEL_MAP

BASE = Path('/mnt/data')
OUT = BASE / 'figure4_molecular_tree_regenerated.png'

# Read tree and prune to the six focal taxa.
tree = Phylo.read(TREE_FILE, 'nexus')
keep_labels = set(TREE_LABEL_MAP.values())
for term in list(tree.get_terminals()):
    if term.name not in keep_labels:
        tree.prune(term)

# Relabel tips to manuscript labels.
reverse_map = {v: k for k, v in TREE_LABEL_MAP.items()}
for term in tree.get_terminals():
    term.name = reverse_map.get(term.name, term.name)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)
Phylo.draw(tree, do_show=False, axes=ax)
ax.set_title('Molecular Phylogenetic Tree of Six Solanaceae Species')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
