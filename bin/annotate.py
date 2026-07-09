#!/usr/bin/env python3
"""
annotate.py — automatic broad cell-type annotation from a marker config.

Given the WNN MuData:
  * Leiden-cluster on the WNN graph at one or more resolutions,
  * score each marker set (sc.tl.score_genes) on the RNA modality,
  * (optionally) also score gene-activity marker sets on the ATAC side,
  * assign each *cluster* the cell type with the highest mean score, using a
    margin rule so ambiguous clusters become "unknown" instead of a coin-flip,
  * write per-cell labels + a cluster->celltype table + a score matrix.

Marker sets come from a YAML you control (assets/markers.yaml), so adding a
lineage or a region set never means touching the code.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import yaml

import mwnn_utils as U


def load_markers(path):
    spec = yaml.safe_load(Path(path).read_text())
    # spec: {celltypes: {name: {rna: [...], atac_genes: [...]}}, params: {...}}
    return spec["celltypes"], spec.get("params", {})


def score_sets(adata, celltypes, key):
    """Add one score column per cell type using genes present in ``adata``."""
    present = set(adata.var_names)
    scored = {}
    for ct, spec in celltypes.items():
        genes = [g for g in spec.get(key, []) if g in present]
        if not genes:
            continue
        col = f"score_{ct}"
        sc.tl.score_genes(adata, genes, score_name=col, use_raw=False)
        scored[ct] = col
    return scored


def assign_clusters(obs, cluster_col, score_cols, margin, min_score):
    """Per-cluster argmax over mean z-scored marker scores."""
    scores = obs[list(score_cols.values())].copy()
    scores.columns = list(score_cols.keys())
    # z-score each marker set across cells so lineages are comparable
    z = (scores - scores.mean()) / (scores.std(ddof=0) + 1e-9)
    z[cluster_col] = obs[cluster_col].values
    cluster_mean = z.groupby(cluster_col, observed=True).mean()

    labels, conf = {}, {}
    for cl, row in cluster_mean.iterrows():
        order = row.sort_values(ascending=False)
        top, second = order.index[0], (order.index[1] if len(order) > 1 else None)
        top_val = order.iloc[0]
        gap = top_val - (order.iloc[1] if len(order) > 1 else -np.inf)
        if top_val < min_score or gap < margin:
            labels[cl] = "unknown"
        else:
            labels[cl] = top
        conf[cl] = {"top": top, "top_score": float(top_val),
                    "second": second, "margin": float(gap)}
    return labels, conf, cluster_mean


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mdata", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--gene_activity", default=None,
                    help="optional ATAC gene-activity h5ad for cross-check")
    ap.add_argument("--resolutions", default="0.05,0.1,0.2,0.5,1.0,1.5,2.0,3.0",)
    ap.add_argument("--primary_resolution", type=float, default=0.2)
    ap.add_argument("--margin", type=float, default=0.25)
    ap.add_argument("--min_score", type=float, default=0.0)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    celltypes, params = load_markers(args.markers)
    margin = params.get("margin", args.margin)
    min_score = params.get("min_score", args.min_score)

    mdata = mu.read(args.mdata)

    # Leiden on the WNN graph at each requested resolution
    res_list = [float(r) for r in args.resolutions.split(",")]
    for r in res_list:
        key = f"leiden_wnn_{str(r).replace('.', '')}"
        sc.tl.leiden(mdata, neighbors_key="neighbors", resolution=r,
                     use_weights=True, directed=False, key_added=key)
    primary_key = f"leiden_wnn_{str(args.primary_resolution).replace('.', '')}"

    # Score marker sets on RNA
    rna = mdata["rna"]
    rna_scores = score_sets(rna, celltypes, key="rna")
    # propagate scores + cluster labels to a single obs frame
    obs = rna.obs.copy()
    obs[primary_key] = mdata.obs[primary_key].values

    labels, conf, cluster_mean = assign_clusters(
        obs, primary_key, rna_scores, margin=margin, min_score=min_score
    )

    mdata.obs["broad_celltype"] = (
        mdata.obs[primary_key].astype(str).map(labels).astype("category")
    )
    # carry scores up for plotting
    for ct, col in rna_scores.items():
        mdata.obs[col] = rna.obs[col].values

    # Optional ATAC gene-activity cross-check (advisory only)
    if args.gene_activity and Path(args.gene_activity).exists():
        ga = sc.read_h5ad(args.gene_activity)
        ga = ga[ga.obs_names.isin(mdata.obs_names)].copy()
        ga = ga[mdata.obs_names.intersection(ga.obs_names)].copy()
        atac_scores = score_sets(ga, celltypes, key="atac_genes")
        if atac_scores:
            ga_obs = ga.obs.copy()
            ga_obs[primary_key] = mdata.obs.loc[ga.obs_names, primary_key].values
            _, atac_conf, atac_cluster_mean = assign_clusters(
                ga_obs, primary_key, atac_scores, margin=0.0, min_score=-np.inf
            )
            atac_cluster_mean.to_csv(out / "atac_cluster_scores.csv")

    # persist tables
    cluster_mean.to_csv(out / "rna_cluster_scores.csv")
    pd.DataFrame([
        {"cluster": cl, "broad_celltype": labels[cl], **conf[cl]}
        for cl in labels
    ]).to_csv(out / "cluster_assignments.csv", index=False)

    ann_cols = ["sample", "broad_celltype", primary_key] + \
               [f"leiden_wnn_{str(r).replace('.', '')}" for r in res_list]
    mdata.obs[ann_cols].to_csv(out / "cell_annotations.tsv", sep="\t",
                               index_label="cell_id")

    mdata.write_h5mu(args.out)
    vc = mdata.obs["broad_celltype"].value_counts()
    U.log("broad_celltype counts:\n" + vc.to_string())


if __name__ == "__main__":
    main()
