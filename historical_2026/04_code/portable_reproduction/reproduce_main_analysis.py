from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[2]
IMAGE_LEVEL_FILE = ROOT / "03_processed_data" / "image_level" / "solanaceae_image_level_medians.csv"
MOLECULAR_DISTANCE_FILE = (
    ROOT
    / "03_processed_data"
    / "distance_matrices"
    / "solanaceae_molecular_distance_matrix_from_tree.csv"
)
OUT_DIR = ROOT / "05_results" / "reproduced_tables"

SPECIES_ORDER = [
    "B. suaveolens",
    "C. annuum var. grossum",
    "L. chinense Mill.",
    "P. × hybrida ‘Dreams Red’",
    "P. alkekengi var. franchetii",
    "S. nigrum",
]

FOCAL_TRAITS = ["tortuosity", "curvature_mean_rad", "curv_per_len"]


def zscore(df: pd.DataFrame) -> pd.DataFrame:
    return (df - df.mean(axis=0)) / df.std(axis=0, ddof=0)


def upper_triangle_values(mat: np.ndarray) -> np.ndarray:
    idx = np.triu_indices_from(mat, k=1)
    return mat[idx]


def mantel_correlation(x_mat: np.ndarray, y_mat: np.ndarray, method: str) -> float:
    x = upper_triangle_values(x_mat)
    y = upper_triangle_values(y_mat)
    if method == "spearman":
        return float(spearmanr(x, y).statistic)
    if method == "pearson":
        return float(pearsonr(x, y).statistic)
    raise ValueError("method must be 'spearman' or 'pearson'")


def exact_mantel(
    x_mat: np.ndarray,
    y_mat: np.ndarray,
    method: str = "spearman",
    alternative: str = "greater",
) -> tuple[float, float]:
    obs = mantel_correlation(x_mat, y_mat, method=method)
    n = x_mat.shape[0]
    perm_stats = []
    for perm in itertools.permutations(range(n)):
        perm = np.asarray(perm)
        y_perm = y_mat[np.ix_(perm, perm)]
        perm_stats.append(mantel_correlation(x_mat, y_perm, method=method))
    perm_stats = np.asarray(perm_stats, dtype=float)
    if alternative == "greater":
        p_value = float(np.mean(perm_stats >= obs - 1e-12))
    elif alternative == "two-sided":
        p_value = float(np.mean(np.abs(perm_stats) >= abs(obs) - 1e-12))
    else:
        raise ValueError("alternative must be 'greater' or 'two-sided'")
    return obs, p_value


def load_image_level(use_contraction: bool = False) -> pd.DataFrame:
    df = pd.read_csv(IMAGE_LEVEL_FILE)
    df["curv_per_len"] = np.log1p(df["curvature_sum_rad"] / df["length_um"])
    trait_cols = FOCAL_TRAITS.copy()
    if use_contraction:
        df["contraction_abs"] = np.log1p(np.abs(df["contraction_ratio"]))
        trait_cols.append("contraction_abs")
    df[trait_cols] = zscore(df[trait_cols].astype(float))
    return df


def species_centroids(image_df: pd.DataFrame, use_contraction: bool = False) -> pd.DataFrame:
    trait_cols = FOCAL_TRAITS.copy()
    if use_contraction:
        trait_cols.append("contraction_abs")
    centroids = image_df.groupby("species")[trait_cols].median(numeric_only=True)
    return centroids.loc[SPECIES_ORDER]


def distance_matrix(centroids: pd.DataFrame) -> pd.DataFrame:
    mat = squareform(pdist(centroids.values, metric="euclidean"))
    return pd.DataFrame(mat, index=centroids.index, columns=centroids.index)


def load_molecular_distance() -> pd.DataFrame:
    mol = pd.read_csv(MOLECULAR_DISTANCE_FILE, index_col=0)
    return mol.loc[SPECIES_ORDER, SPECIES_ORDER]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    image_df = load_image_level(use_contraction=False)
    centroids = species_centroids(image_df, use_contraction=False)
    morph = distance_matrix(centroids)
    mol = load_molecular_distance()

    rows = []
    for method in ["spearman", "pearson"]:
        r, p = exact_mantel(morph.values, mol.values, method=method, alternative="greater")
        rows.append({"analysis": "main_focal_block", "method": method, "r_M": r, "p_exact": p})

    image_df_contraction = load_image_level(use_contraction=True)
    centroids_contraction = species_centroids(image_df_contraction, use_contraction=True)
    morph_contraction = distance_matrix(centroids_contraction)
    for method in ["spearman", "pearson"]:
        r, p = exact_mantel(
            morph_contraction.values,
            mol.values,
            method=method,
            alternative="greater",
        )
        rows.append({"analysis": "focal_plus_contraction", "method": method, "r_M": r, "p_exact": p})

    loo_rows = []
    for excluded in SPECIES_ORDER:
        kept = [sp for sp in SPECIES_ORDER if sp != excluded]
        morph_sub = morph.loc[kept, kept]
        mol_sub = mol.loc[kept, kept]
        r, p = exact_mantel(morph_sub.values, mol_sub.values, method="spearman", alternative="greater")
        loo_rows.append({"excluded_species": excluded, "method": "spearman", "r_M": r, "p_exact": p})

    centroids.to_csv(OUT_DIR / "reproduced_species_centroids_path_curvature.csv")
    morph.to_csv(OUT_DIR / "reproduced_morph_distance_matrix_path_curvature.csv")
    pd.DataFrame(rows).to_csv(OUT_DIR / "reproduced_mantel_summary.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(OUT_DIR / "reproduced_leave_one_species_out.csv", index=False)

    print("Main focal block")
    print(pd.DataFrame(rows).to_string(index=False))
    print()
    print("Leave-one-species-out Spearman exact Mantel")
    print(pd.DataFrame(loo_rows).to_string(index=False))


if __name__ == "__main__":
    main()
