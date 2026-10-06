#!/usr/bin/env python3
"""
Xen-Spat: size of the 25-cell niche neighborhood in micrometers.

The niche analysis (Seurat BuildNicheAssay, neighbors.k = 25) defines each
cell's neighborhood as its 25 nearest cells, counting the cell itself, as
Seurat's k-nearest-neighbor search does. This script reports how large that
neighborhood is in micrometers. For each cell it measures the Euclidean
distance from the cell's centroid to the farthest cell in the neighborhood,
which is the 24th nearest other cell. It then summarizes these distances per
sample as the median and interquartile range.

Input (one of the following):
  --h5ad  AnnData file with obs columns x_centroid, y_centroid (um) and batch
          (e.g. Integrated_Data_CelltypeAnnotated_9225.h5ad). Batch codes are
          mapped to sample IDs with BATCH_TO_SAMPLE below.
  --csv   CSV with columns sample, x_centroid, y_centroid (um), one row per
          cell, e.g. the uniform label files concatenated across samples.

Cells are those retained after quality filtering (Table S3). Coordinates are
the Xenium cell centroids in um.

Output: a CSV with one row per sample (n_cells, median, Q1, Q3 in um), also
printed to the screen.

Usage:
  python niche_neighborhood_radius.py --h5ad Integrated_Data_CelltypeAnnotated_9225.h5ad \
      --out niche_neighborhood_radius_by_sample.csv

Requirements: numpy, pandas, scipy (and anndata if --h5ad is used).
"""
import argparse

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

K_NEIGHBORS = 25  # neighbors.k used in BuildNicheAssay (includes the cell itself)

BATCH_TO_SAMPLE = {"BC1": "046", "BC2": "047", "BC3": "049",
                   "BC4": "089", "BC5": "097", "BC6": "100"}


def load_cells(h5ad=None, csv=None):
    if h5ad:
        import anndata as ad
        obs = ad.read_h5ad(h5ad, backed="r").obs
        df = obs[["x_centroid", "y_centroid", "batch"]].copy()
        df["sample"] = df["batch"].astype(str).map(BATCH_TO_SAMPLE)
        if df["sample"].isna().any():
            raise ValueError("Unrecognized batch codes: "
                             f"{sorted(df.loc[df['sample'].isna(), 'batch'].unique())}")
    else:
        df = pd.read_csv(csv, dtype={"sample": str})
        df["sample"] = df["sample"].str.zfill(3) if df["sample"].str.isdigit().all() else df["sample"]
    return df[["sample", "x_centroid", "y_centroid"]]


def neighborhood_radius(xy, k=K_NEIGHBORS):
    """Distance (um) from each cell to the farthest member of its k-cell
    neighborhood, where the neighborhood includes the cell itself."""
    if len(xy) < k:
        raise ValueError(f"Sample has fewer than {k} cells")
    dist, _ = cKDTree(xy).query(xy, k=k)  # column 0 is the cell itself (distance 0)
    return dist[:, k - 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--h5ad")
    src.add_argument("--csv")
    ap.add_argument("--out", default="niche_neighborhood_radius_by_sample.csv")
    ap.add_argument("--k", type=int, default=K_NEIGHBORS)
    a = ap.parse_args()

    cells = load_cells(a.h5ad, a.csv)
    rows = []
    for sample, g in cells.groupby("sample", sort=True):
        r = neighborhood_radius(g[["x_centroid", "y_centroid"]].to_numpy(float), a.k)
        q1, med, q3 = np.percentile(r, [25, 50, 75])
        rows.append({"sample": sample, "n_cells": len(g),
                     "median_radius_um": round(med, 1),
                     "q1_um": round(q1, 1), "q3_um": round(q3, 1)})
    res = pd.DataFrame(rows)
    res.to_csv(a.out, index=False)
    print(res.to_string(index=False))
    print(f"\nRange of per-sample medians: {res.median_radius_um.min()} to "
          f"{res.median_radius_um.max()} um (k = {a.k}, including the cell itself)")


if __name__ == "__main__":
    main()
