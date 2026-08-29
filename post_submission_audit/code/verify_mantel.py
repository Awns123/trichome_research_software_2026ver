"""Recompute exact Mantel statistics from archived and audited matrices.

This is a 2026-08-26 public-release verifier, not historical student code and
not a phylogenetic-tree inference pipeline. It uses all 6! simultaneous taxon-
label permutations, including the identity permutation.
"""

from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

MATRIX_FILES = {
    "archived_tree": "molecular_matrix_original_archived_tree.csv",
    "raw_as_stored": "molecular_matrix_original_as_stored_16taxa_focal6.csv",
    "feature_oriented": (
        "molecular_matrix_feature_oriented_original_accessions_16taxa_focal6.csv"
    ),
    "minimal_corrected": "molecular_matrix_minimal_corrected_17taxa_focal6.csv",
    "standardized_hybrid": "molecular_matrix_standardized_hybrid_17taxa_focal6.csv",
}

RESULT_COLUMNS = [
    "analysis",
    "taxa",
    "exact_label_permutations",
    "spearman_r",
    "spearman_p_greater",
    "spearman_p_two_sided",
    "pearson_r",
    "pearson_p_greater",
    "pearson_p_two_sided",
]


def load_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, index_col=0)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Matrix is not square: {path}")
    if list(matrix.index) != list(matrix.columns):
        raise ValueError(f"Row and column labels differ: {path}")
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite value found: {path}")
    if not np.allclose(values, values.T, atol=1e-12, rtol=0):
        raise ValueError(f"Matrix is not symmetric: {path}")
    if not np.allclose(np.diag(values), 0, atol=1e-12, rtol=0):
        raise ValueError(f"Diagonal is not zero: {path}")
    return matrix


def upper_triangle(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices_from(matrix, k=1)]


def correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    x_upper = upper_triangle(x)
    y_upper = upper_triangle(y)
    if method == "spearman":
        return float(spearmanr(x_upper, y_upper).statistic)
    if method == "pearson":
        return float(pearsonr(x_upper, y_upper).statistic)
    raise ValueError(f"Unsupported method: {method}")


def exact_test(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float, float]:
    observed = correlation(x, y, method)
    permuted = []
    for permutation in itertools.permutations(range(x.shape[0])):
        order = np.asarray(permutation)
        y_permuted = y[np.ix_(order, order)]
        permuted.append(correlation(x, y_permuted, method))
    permuted_array = np.asarray(permuted, dtype=float)
    tolerance = 1e-12
    p_greater = float(np.mean(permuted_array >= observed - tolerance))
    p_two_sided = float(
        np.mean(np.abs(permuted_array) >= abs(observed) - tolerance)
    )
    return observed, p_greater, p_two_sided


def recompute() -> pd.DataFrame:
    morphology = load_matrix(DATA / "morphology_matrix_original.csv")
    rows: list[dict[str, float | int | str]] = []
    for analysis, filename in MATRIX_FILES.items():
        molecular = load_matrix(DATA / filename)
        if set(morphology.index) != set(molecular.index):
            raise ValueError(f"Taxon labels differ for {analysis}")
        molecular = molecular.loc[morphology.index, morphology.index]
        spearman = exact_test(
            morphology.to_numpy(dtype=float),
            molecular.to_numpy(dtype=float),
            "spearman",
        )
        pearson = exact_test(
            morphology.to_numpy(dtype=float),
            molecular.to_numpy(dtype=float),
            "pearson",
        )
        rows.append(
            {
                "analysis": analysis,
                "taxa": len(morphology),
                "exact_label_permutations": math.factorial(len(morphology)),
                "spearman_r": spearman[0],
                "spearman_p_greater": spearman[1],
                "spearman_p_two_sided": spearman[2],
                "pearson_r": pearson[0],
                "pearson_p_greater": pearson[1],
                "pearson_p_two_sided": pearson[2],
            }
        )
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


def compare_with_archived(recomputed: pd.DataFrame, tolerance: float = 1e-12) -> None:
    archived = pd.read_csv(DATA / "mantel_all_sensitivity.csv").set_index("analysis")
    recomputed_indexed = recomputed.set_index("analysis")
    numeric_columns = [
        "spearman_r",
        "spearman_p_greater",
        "spearman_p_two_sided",
        "pearson_r",
        "pearson_p_greater",
        "pearson_p_two_sided",
    ]
    missing = sorted(set(recomputed_indexed.index) - set(archived.index))
    if missing:
        raise AssertionError(f"Archived summary is missing rows: {missing}")
    for analysis in recomputed_indexed.index:
        for column in numeric_columns:
            observed = float(recomputed_indexed.loc[analysis, column])
            expected = float(archived.loc[analysis, column])
            if not np.isclose(observed, expected, atol=tolerance, rtol=0):
                raise AssertionError(
                    f"Mismatch for {analysis}.{column}: {observed} != {expected}"
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if recomputed values differ from the archived audit summary.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional CSV output path. Parent directories are created.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = recompute()
    if args.check:
        compare_with_archived(result)
        print("CHECK_OK: recomputed statistics match the archived audit summary")
    print(result.to_string(index=False))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
