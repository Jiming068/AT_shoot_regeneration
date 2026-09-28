#!/usr/bin/env python3
"""Classify SIM6 sub-1 niches and surrounding sub-3/cluster-10 domains.

The implementation follows the final analysis used for the manuscript figure:
1. Detect valid sub-1 niches with DBSCAN (eps=50, min_samples=5).
2. Within each niche, split sub-1 cells by 1D K-means on distance to the
   coordinate-wise median centroid; the nearer partition is core.
3. For sub-3 and cluster-10 cells, calculate distance to the nearest valid
   sub-1 cell. Distances <=200 are shell; distances >200 are displaced.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN, KMeans

TARGET_TIME = "SIM6"
REFERENCE = "sub-1"
RING_TYPES = ("sub-3", "10")
EPS = 50.0
MIN_SAMPLES = 5
SHELL_RADIUS = 200.0
RANDOM_SEED = 123

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Metadata CSV or CSV.GZ from step 01")
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.input, dtype={"cell_id": str, "time": str, "cluster": str, "tissue": str})
    required = {"cell_id", "time", "cluster", "tissue", "x", "y"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if meta["cell_id"].duplicated().any():
        raise ValueError("cell_id must be unique")

    target = meta.loc[
        (meta["time"] == TARGET_TIME) & meta["cluster"].isin((REFERENCE,) + RING_TYPES)
    ].copy()
    target["niche_id"] = pd.array([pd.NA] * len(target), dtype="Int64")
    target["distance_to_niche_centroid"] = np.nan
    target["distance_to_nearest_valid_sub1"] = np.nan
    target["niche_class"] = pd.NA

    ref_mask = target["cluster"].eq(REFERENCE).to_numpy()
    ref_coords = target.loc[ref_mask, ["x", "y"]].to_numpy(float)
    model = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES)
    labels = model.fit_predict(ref_coords)
    target.loc[ref_mask, "niche_id"] = pd.array(
        [value if value >= 0 else pd.NA for value in labels], dtype="Int64"
    )

    ref_positions = np.flatnonzero(ref_mask)
    for niche_id in sorted(set(labels).difference({-1})):
        local = np.flatnonzero(labels == niche_id)
        positions = ref_positions[local]
        coords = ref_coords[local]
        centroid = np.median(coords, axis=0)
        distances = np.linalg.norm(coords - centroid, axis=1)
        km = KMeans(n_clusters=2, random_state=RANDOM_SEED, n_init=10)
        raw_group = km.fit_predict(distances.reshape(-1, 1))
        mean_distance = {
            group: float(distances[raw_group == group].mean()) for group in np.unique(raw_group)
        }
        core_group = min(mean_distance, key=mean_distance.get)
        target.iloc[positions, target.columns.get_loc("distance_to_niche_centroid")] = distances
        target.iloc[positions, target.columns.get_loc("niche_class")] = np.where(
            raw_group == core_group, "sub-1 core", "sub-1 bdry"
        )

    valid_ref_coords = ref_coords[labels >= 0]
    if valid_ref_coords.size == 0:
        raise RuntimeError("DBSCAN detected no valid sub-1 niche")
    nonref_positions = np.flatnonzero(~ref_mask)
    nonref_coords = target.iloc[nonref_positions][["x", "y"]].to_numpy(float)
    nearest_distance = cKDTree(valid_ref_coords).query(nonref_coords, k=1)[0]
    target.iloc[
        nonref_positions, target.columns.get_loc("distance_to_nearest_valid_sub1")
    ] = nearest_distance

    nonref_types = target.iloc[nonref_positions]["cluster"].to_numpy(str)
    nonref_class = np.where(
        nonref_types == "sub-3",
        np.where(nearest_distance <= SHELL_RADIUS, "sub-3 shell", "sub-3 disp."),
        np.where(nearest_distance <= SHELL_RADIUS, "10 shell", "10 disp."),
    )
    target.iloc[nonref_positions, target.columns.get_loc("niche_class")] = nonref_class

    output_columns = [
        "cell_id", "time", "cluster", "tissue", "x", "y", "niche_id",
        "distance_to_niche_centroid", "distance_to_nearest_valid_sub1", "niche_class",
    ]
    classification_file = outdir / "SIM6_sub1_sub3_10_niche_classification_reproduced.csv"
    target[output_columns].to_csv(classification_file, index=False)

    summary = (
        target["niche_class"].fillna("unclassified sub-1 noise").value_counts(dropna=False)
        .rename_axis("niche_class").reset_index(name="n_cells")
    )
    summary.to_csv(outdir / "SIM6_niche_class_summary_reproduced.csv", index=False)

    parameters = {
        "target_time": TARGET_TIME,
        "reference_celltype": REFERENCE,
        "ring_celltypes": list(RING_TYPES),
        "dbscan_eps": EPS,
        "dbscan_min_samples": MIN_SAMPLES,
        "shell_radius": SHELL_RADIUS,
        "kmeans_random_seed": RANDOM_SEED,
        "kmeans_n_init": 10,
        "n_target_cells": int(len(target)),
        "n_valid_niches": int(len(set(labels).difference({-1}))),
        "n_unclassified_sub1_noise": int((labels == -1).sum()),
        "class_counts": dict(zip(summary["niche_class"], summary["n_cells"].astype(int))),
    }
    with open(outdir / "SIM6_niche_reclassification_parameters.json", "w") as stream:
        json.dump(parameters, stream, indent=2)

    print(summary.to_string(index=False))
    print(f"Saved: {classification_file}")

if __name__ == "__main__":
    main()
