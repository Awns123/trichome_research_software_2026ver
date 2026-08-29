
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway

BASE = Path("/mnt/data")
IMAGE_LEVEL_FILE = BASE / "solanaceae_image_level_medians.csv"

def permanova_from_dm(dm: np.ndarray, groups, n_perm: int = 10000, seed: int = 123):
    dm = np.asarray(dm, dtype=float)
    groups = np.asarray(groups)
    n = dm.shape[0]
    uniq = pd.unique(groups)
    g = len(uniq)

    A = -0.5 * (dm ** 2)
    H = np.eye(n) - np.ones((n, n)) / n
    G = H @ A @ H
    ss_total = np.trace(G)

    def ss_within(groups_local):
        ssw = 0.0
        for gr in pd.unique(groups_local):
            idx = np.where(groups_local == gr)[0]
            ng = len(idx)
            if ng <= 1:
                continue
            Ag = A[np.ix_(idx, idx)]
            Hg = np.eye(ng) - np.ones((ng, ng)) / ng
            Gg = Hg @ Ag @ Hg
            ssw += np.trace(Gg)
        return float(ssw)

    ssw = ss_within(groups)
    ssb = ss_total - ssw
    dfb = g - 1
    dfw = n - g
    F_obs = (ssb / dfb) / (ssw / dfw)

    rng = np.random.default_rng(seed)
    perm_F = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        gp = rng.permutation(groups)
        ssw_p = ss_within(gp)
        ssb_p = ss_total - ssw_p
        perm_F[i] = (ssb_p / dfb) / (ssw_p / dfw)

    p_value = (np.sum(perm_F >= F_obs) + 1) / (n_perm + 1)

    return {
        "F": float(F_obs),
        "p_value": float(p_value),
        "R2": float(ssb / ss_total),
        "df_between": int(dfb),
        "df_within": int(dfw),
        "ss_between": float(ssb),
        "ss_within": float(ssw),
    }

def pcoa_coords(dm: np.ndarray, tol: float = 1e-12):
    dm = np.asarray(dm, dtype=float)
    n = dm.shape[0]
    A = -0.5 * (dm ** 2)
    H = np.eye(n) - np.ones((n, n)) / n
    G = H @ A @ H
    eigvals, eigvecs = np.linalg.eigh(G)
    order = eigvals.argsort()[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = eigvals > tol
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    return coords

def permdisp_from_dm(dm: np.ndarray, groups, n_perm: int = 10000, seed: int = 123):
    dm = np.asarray(dm, dtype=float)
    groups = np.asarray(groups)

    coords = pcoa_coords(dm)
    uniq = pd.unique(groups)

    dists = np.empty(len(groups), dtype=float)
    for gr in uniq:
        idx = np.where(groups == gr)[0]
        centroid = coords[idx].mean(axis=0, keepdims=True)
        dists[idx] = np.sqrt(((coords[idx] - centroid) ** 2).sum(axis=1))

    obs_F = float(f_oneway(*[dists[groups == gr] for gr in uniq]).statistic)

    rng = np.random.default_rng(seed)
    perm_F = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        gp = rng.permutation(groups)
        perm_F[i] = float(f_oneway(*[dists[gp == gr] for gr in uniq]).statistic)

    p_value = (np.sum(perm_F >= obs_F) + 1) / (n_perm + 1)

    return {
        "F": obs_F,
        "p_value": float(p_value),
        "group_mean_distances": {
            str(gr): float(dists[groups == gr].mean()) for gr in uniq
        },
    }

if __name__ == "__main__":
    df = pd.read_csv(IMAGE_LEVEL_FILE).copy()

    # focal trait block used in the manuscript
    df["curv_per_len"] = df["curvature_sum_rad"] / df["length_um"]
    df["curv_per_len"] = np.log1p(df["curv_per_len"])

    traits = ["tortuosity", "curvature_mean_rad", "curv_per_len"]
    X = df[traits].astype(float).copy()
    X = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)

    groups = df["species"].astype(str).values
    dm = squareform(pdist(X.values, metric="euclidean"))

    permanova_res = permanova_from_dm(dm, groups, n_perm=10000, seed=123)
    permdisp_res = permdisp_from_dm(dm, groups, n_perm=10000, seed=123)

    print("=== PERMANOVA ===")
    for k, v in permanova_res.items():
        print(f"{k}: {v}")

    print("\n=== PERMDISP ===")
    print(f"F: {permdisp_res['F']}")
    print(f"p_value: {permdisp_res['p_value']}")
    print("group_mean_distances:")
    for k, v in permdisp_res["group_mean_distances"].items():
        print(f"  {k}: {v}")
