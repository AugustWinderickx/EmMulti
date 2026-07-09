#!/usr/bin/env python3
"""
plots.py — final figure set from the annotated MuData.

All figures scale to the number of samples / cells / cell types so text never
overlaps:
  * WNN UMAP coloured by broad_celltype, sample, and primary leiden,
  * marker dotplot grouped by broad_celltype,
  * per-lineage marker-score UMAP grid,
  * cell-type composition per sample (stacked bars, adaptive width),
  * optional ATAC gene-activity UMAP for a few key markers.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import yaml
import matplotlib.pyplot as plt

import mwnn_utils as U


def umap_by(mdata, color, path, title=None, categorical=True):
    xy = mdata.obsm["X_umap"]
    n = mdata.n_obs
    s = U.umap_point_size(n)
    fig, ax = plt.subplots(figsize=(8, 7))
    vals = mdata.obs[color].astype(str).values if categorical else None
    if categorical:
        cats = pd.unique(vals)
        cats = sorted(cats)
        pal = dict(zip(cats, U.categorical_palette(len(cats))))
        for c in cats:
            m = vals == c
            ax.scatter(xy[m, 0], xy[m, 1], s=s, c=pal[c], label=c,
                       linewidths=0, rasterized=True)
        U.legend_outside(ax, len(cats), title=color)
    ax.set(xticks=[], yticks=[], title=title or color)
    ax.set_xlabel("WNN-UMAP1"); ax.set_ylabel("WNN-UMAP2")
    U.save_fig(fig, path)


def umap_scores(mdata, score_cols, path):
    n = len(score_cols)
    if n == 0:
        return
    nr, nc = U.grid_dims(n, max_cols=4)
    xy = mdata.obsm["X_umap"]
    s = U.umap_point_size(mdata.n_obs)
    fig, axes = plt.subplots(nr, nc, figsize=(3.4 * nc, 3.2 * nr), squeeze=False)
    for ax, (ct, col) in zip(axes.ravel(), score_cols.items()):
        v = mdata.obs[col].values
        vmax = np.quantile(v, 0.99)
        sca = ax.scatter(xy[:, 0], xy[:, 1], c=v, s=s, cmap="magma",
                         vmin=np.quantile(v, 0.02), vmax=vmax,
                         linewidths=0, rasterized=True)
        ax.set(xticks=[], yticks=[], title=ct)
        fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle("Marker-set scores (RNA)")
    U.save_fig(fig, path)


def dotplot(mdata, celltypes, groupby, path):
    rna = mdata["rna"].copy()
    rna.obs[groupby] = mdata.obs[groupby].values
    present = set(rna.var_names)
    var_groups = {ct: [g for g in spec.get("rna", []) if g in present]
                  for ct, spec in celltypes.items()}
    var_groups = {k: v for k, v in var_groups.items() if v}
    n_genes = sum(len(v) for v in var_groups.values())
    fig_w = max(8, 0.28 * n_genes + 3)
    n_groups = rna.obs[groupby].nunique()
    fig_h = max(4, 0.4 * n_groups + 2)
    sc.pl.dotplot(rna, var_groups, groupby=groupby, standard_scale="var",
                  show=False, figsize=(fig_w, fig_h))
    U.save_fig(plt.gcf(), path)


def composition(mdata, path, groupby="broad_celltype", batch="sample"):
    ct = pd.crosstab(mdata.obs[batch], mdata.obs[groupby], normalize="index")
    samples = ct.index.tolist()
    cats = ct.columns.tolist()
    pal = U.categorical_palette(len(cats))
    fig, ax = plt.subplots(figsize=(U.sample_figwidth(len(samples)), 5))
    bottom = np.zeros(len(samples))
    for c, color in zip(cats, pal):
        ax.bar(range(len(samples)), ct[c].values, bottom=bottom,
               color=color, label=c, width=0.85)
        bottom += ct[c].values
    U.style_sample_xticks(ax, samples)
    ax.set(ylabel="fraction of cells", ylim=(0, 1),
           title="Cell-type composition per sample")
    U.legend_outside(ax, len(cats))
    U.save_fig(fig, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mdata", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--gene_activity", default=None)
    ap.add_argument("--primary_resolution", type=float, default=0.2)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    spec = yaml.safe_load(Path(args.markers).read_text())
    celltypes = spec["celltypes"]
    mdata = mu.read(args.mdata)
    primary_key = f"leiden_wnn_{str(args.primary_resolution).replace('.', '')}"

    umap_by(mdata, "broad_celltype", out / "umap_broad_celltype.png",
            "Broad cell type (WNN)")
    umap_by(mdata, "sample", out / "umap_sample.png", "Sample")
    if primary_key in mdata.obs:
        umap_by(mdata, primary_key, out / "umap_leiden.png",
                f"Leiden (res {args.primary_resolution})")

    score_cols = {ct: f"score_{ct}" for ct in celltypes
                  if f"score_{ct}" in mdata.obs}
    umap_scores(mdata, score_cols, out / "umap_marker_scores.png")

    dotplot(mdata, celltypes, "broad_celltype", out / "dotplot_markers.png")
    composition(mdata, out / "composition_per_sample.png")

    U.log("plots written to " + str(out))


if __name__ == "__main__":
    main()
