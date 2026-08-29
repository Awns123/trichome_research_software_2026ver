
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import shapiro, kruskal

BASE = Path("/mnt/data")

SPECIES_FILES = {
    "B. suaveolens": BASE / "B. suaveolens.csv",
    "C. annuum var. grossum": BASE / "C. annuum var. grossum.csv",
    "L. chinense Mill.": BASE / "L. chinense Mill..csv",
    "P. × hybrida ‘Dreams Red’": BASE / "P. × hybrida ‘Dreams Red’.csv",
    "P. alkekengi var. franchetii": BASE / "P. alkekengi var. franchetii.csv",
    "S. nigrum": BASE / "S. nigrum.csv",
}

TRAIT_COLS = [
    "area_um2",
    "perimeter_um",
    "bbox_aspect_ratio",
    "ellipse_major_um",
    "ellipse_minor_um",
    "ellipse_angle_deg",
    "contraction_ratio",
    "num_endpoints",
    "num_branchpoints",
    "length_um",
    "straight_length_um",
    "tortuosity",
    "width_base_um",
    "width_tip_um",
    "width_max_um",
    "curvature_mean_rad",
    "curvature_max_rad",
    "curvature_sum_rad",
    "skeleton_length_um",
    "skeleton_ratio",
    "upper_curv_mean_rad",
    "upper_curv_max_rad",
    "upper_curv_sum_rad",
    "lower_curv_mean_rad",
    "lower_curv_max_rad",
    "lower_curv_sum_rad",
]

def build_image_level_medians() -> pd.DataFrame:
    tables = []

    for species, path in SPECIES_FILES.items():
        df = pd.read_csv(path)

        # image-level median
        img = (
            df.groupby("image", dropna=False)[TRAIT_COLS]
            .median(numeric_only=True)
            .reset_index()
        )

        # derived variable used later in the manuscript
        img["curv_per_len"] = img["curvature_sum_rad"] / img["length_um"]
        img["species"] = species
        tables.append(img)

    all_img = pd.concat(tables, ignore_index=True)
    return all_img


def shapiro_wilk_by_trait(image_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for trait in TRAIT_COLS:
        x = image_df[trait].dropna().astype(float).values
        W, p = shapiro(x)
        rows.append({
            "trait": trait,
            "n": len(x),
            "W": W,
            "p_value": p,
            "normal_at_0_05": bool(p >= 0.05),
        })

    # derived variable tested separately
    x = image_df["curv_per_len"].dropna().astype(float).values
    W, p = shapiro(x)
    rows.append({
        "trait": "curv_per_len",
        "n": len(x),
        "W": W,
        "p_value": p,
        "normal_at_0_05": bool(p >= 0.05),
    })

    return pd.DataFrame(rows).sort_values("p_value")


def kruskal_wallis_by_trait(image_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    species_order = list(SPECIES_FILES.keys())

    for trait in TRAIT_COLS + ["curv_per_len"]:
        groups = []
        for sp in species_order:
            vals = (
                image_df.loc[image_df["species"] == sp, trait]
                .dropna()
                .astype(float)
                .values
            )
            groups.append(vals)

        H, p = kruskal(*groups)
        rows.append({
            "trait": trait,
            "H": H,
            "p_value": p,
            "significant_at_0_05": bool(p < 0.05),
        })

    out = pd.DataFrame(rows).sort_values("p_value")
    return out


if __name__ == "__main__":
    image_df = build_image_level_medians()

    shapiro_df = shapiro_wilk_by_trait(image_df)
    kw_df = kruskal_wallis_by_trait(image_df)

    print("=== image-level median dataset ===")
    print(f"n_images = {len(image_df)}")
    print()

    print("=== Shapiro–Wilk results ===")
    print(shapiro_df.to_string(index=False))
    print()

    print("=== Kruskal–Wallis results ===")
    print(kw_df.to_string(index=False))
    print()

    print("=== focal traits ===")
    focal = kw_df[kw_df["trait"].isin(["tortuosity", "curvature_mean_rad", "curv_per_len"])]
    print(focal.to_string(index=False))

    shapiro_df.to_csv(BASE / "solanaceae_image_level_shapiro_wilk_results.csv", index=False)
    kw_df.to_csv(BASE / "solanaceae_image_level_kruskal_wallis_results.csv", index=False)
    image_df.to_csv(BASE / "solanaceae_image_level_medians.csv", index=False)
