from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from reproduce_main_analysis import (
    FOCAL_TRAITS,
    OUT_DIR,
    SPECIES_ORDER,
    distance_matrix,
    load_image_level,
    load_molecular_distance,
    species_centroids,
)


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "05_results" / "reproduced_figures"


def make_pca_biplot(centroids: pd.DataFrame) -> None:
    pca = PCA(n_components=2)
    coords = pca.fit_transform(centroids[FOCAL_TRAITS].values)
    loadings = pca.components_.T

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(coords[:, 0], coords[:, 1])
    for i, species in enumerate(centroids.index):
        ax.text(coords[i, 0], coords[i, 1], species, fontsize=8)

    scale = np.max(np.abs(coords)) * 0.75
    labels = ["Tortuosity", "Mean Curvature", "Curvature per Length"]
    for i, label in enumerate(labels):
        ax.arrow(
            0,
            0,
            loadings[i, 0] * scale,
            loadings[i, 1] * scale,
            head_width=0.04 * scale,
            length_includes_head=True,
        )
        ax.text(loadings[i, 0] * scale * 1.08, loadings[i, 1] * scale * 1.08, label, fontsize=8)

    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}%)")
    ax.set_title("PCA biplot of species-level centroids")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "reproduced_pca_biplot.png", dpi=300)
    plt.close(fig)


def make_distance_heatmaps(morph: pd.DataFrame, mol: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for ax, mat, title in [
        (axes[0], morph, "Morphological distance"),
        (axes[1], mol, "Molecular distance"),
    ]:
        im = ax.imshow(mat.values)
        ax.set_xticks(range(len(mat.columns)))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(mat.index)))
        ax.set_yticklabels(mat.index, fontsize=7)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(FIG_DIR / "reproduced_distance_heatmaps.png", dpi=300)
    plt.close(fig)


def make_leave_one_out_bar() -> None:
    loo_file = OUT_DIR / "reproduced_leave_one_species_out.csv"
    if not loo_file.exists():
        raise FileNotFoundError("Run reproduce_main_analysis.py before reproduce_figures.py")
    loo = pd.read_csv(loo_file)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.barh(loo["excluded_species"], loo["r_M"])
    ax.axvline(0.775, linestyle="--", linewidth=1)
    ax.set_xlabel("Spearman exact Mantel r_M")
    ax.set_title("Leave-one-species-out robustness")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "reproduced_leave_one_species_out.png", dpi=300)
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    image_df = load_image_level(use_contraction=False)
    centroids = species_centroids(image_df, use_contraction=False)
    morph = distance_matrix(centroids)
    mol = load_molecular_distance().loc[SPECIES_ORDER, SPECIES_ORDER]

    make_pca_biplot(centroids)
    make_distance_heatmaps(morph, mol)
    make_leave_one_out_bar()
    print(f"Saved reproduced figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
