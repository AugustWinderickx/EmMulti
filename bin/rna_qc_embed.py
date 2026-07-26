#!/usr/bin/env python3
"""
rna_qc_embed.py — concatenate per-sample RNA, QC, normalise, and embed.

Steps mirror a standard scanpy pipeline but batch-aware and with adaptive QC
figures:
  * concat all samples (outer join),
  * QC metrics (mt / ribo / hb), filtering (--qc_mode basic | strict):
      - basic: min_genes/min_cells (+ optional hard max_pct_mt) only,
      - strict: basic, plus MAD-based outlier removal on counts/genes/
        top20-gene-fraction and a mito MAD + hard-cap filter,
  * scrublet per sample,
  * counts layer, normalize_total + log1p,
  * HVG (batch_key=sample), PCA,
  * Harmony over ``sample`` (shape-safe), neighbours, UMAP.

Outputs the processed RNA .h5ad (with X_pca_harmony) for WNN.
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import stats

import mwnn_utils as U


def mad_outlier(adata, metric, n_mads):
    """Boolean mask of cells more than n_mads MADs from the metric's median."""
    m = adata.obs[metric]
    med = np.median(m)
    mad = stats.median_abs_deviation(m)
    return (m < med - n_mads * mad) | (m > med + n_mads * mad)


def qc_plots(adata, outdir):
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                 jitter=0.4, multi_panel=True, show=False)
    U.save_fig(plt.gcf(), outdir / "rna_qc_violin.png")

    fig, ax = plt.subplots(figsize=(5.5, 5))
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",
                  color="pct_counts_mt", ax=ax, show=False)
    U.save_fig(fig, outdir / "rna_qc_scatter.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad_glob", required=True)
    ap.add_argument("--species", default="human", choices=["human", "mouse"])
    ap.add_argument("--min_genes", type=int, default=100)
    ap.add_argument("--min_cells", type=int, default=3)
    ap.add_argument("--max_pct_mt", type=float, default=None,
                    help="optional hard cap on %% mito; default keep all")
    ap.add_argument("--qc_mode", default="basic", choices=["basic", "strict"],
                    help="basic: min_genes/min_cells (+ optional max_pct_mt) "
                         "only. strict: also drops MAD-based outliers on "
                         "counts/genes/top20-gene-fraction, plus a mito "
                         "MAD + hard-cap filter (sc-best-practices recipe)")
    ap.add_argument("--strict_n_mads", type=float, default=5.0,
                    help="strict mode: MAD threshold for counts/genes/top20")
    ap.add_argument("--strict_mt_n_mads", type=float, default=3.0,
                    help="strict mode: MAD threshold for pct_counts_mt")
    ap.add_argument("--strict_max_pct_mt", type=float, default=8.0,
                    help="strict mode: hard cap on %% mito")
    ap.add_argument("--n_top_genes", type=int, default=3000)
    ap.add_argument("--n_pcs", type=int, default=50)
    ap.add_argument("--harmony_key", default="sample")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(args.h5ad_glob))
    if not files:
        raise SystemExit(f"no RNA h5ads matched {args.h5ad_glob}")

    loaded = []
    for f in files:
        a = sc.read_h5ad(f)
        U.assert_raw_counts(a, f)
        loaded.append(a)
    adata = ad.concat(loaded, join="outer", merge="first", index_unique=None)
    adata.obs_names_make_unique()
    U.log(f"concatenated {adata.n_obs} cells x {adata.n_vars} genes "
          f"from {len(files)} samples")

    mt = "MT-" if args.species == "human" else "mt-"
    ribo = ("RPS", "RPL") if args.species == "human" else ("Rps", "Rpl")
    adata.var["mt"] = adata.var_names.str.startswith(mt)
    adata.var["ribo"] = adata.var_names.str.startswith(ribo)
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]" if args.species == "human"
                                                   else r"^Hb[^(p)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"],
                               percent_top=(20, 50, 100, 200, 500),
                               inplace=True, log1p=True)
    qc_plots(adata, out)

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    if args.max_pct_mt is not None:
        adata = adata[adata.obs["pct_counts_mt"] <= args.max_pct_mt].copy()

    if args.qc_mode == "strict":
        outlier = (
            mad_outlier(adata, "log1p_total_counts", args.strict_n_mads)
            | mad_outlier(adata, "log1p_n_genes_by_counts", args.strict_n_mads)
            | mad_outlier(adata, "pct_counts_in_top_20_genes", args.strict_n_mads)
        )
        mt_outlier = (
            mad_outlier(adata, "pct_counts_mt", args.strict_mt_n_mads)
            | (adata.obs["pct_counts_mt"] > args.strict_max_pct_mt)
        )
        n0 = adata.n_obs
        adata = adata[~(outlier | mt_outlier)].copy()
        U.log(f"strict QC removed {n0 - adata.n_obs} outlier cells "
              f"({int(outlier.sum())} count/gene/top20 outliers, "
              f"{int(mt_outlier.sum())} mito outliers)")

    sc.pp.scrublet(adata, batch_key="sample")
    if "predicted_doublet" in adata.obs:
        n0 = adata.n_obs
        adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
        U.log(f"scrublet removed {n0 - adata.n_obs} doublets")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes,
                                batch_key="sample")
    sc.tl.pca(adata, n_comps=args.n_pcs)

    meta = adata.obs[[args.harmony_key]].copy()
    adata.obsm["X_pca_harmony"] = U.harmony_embedding(
        adata.obsm["X_pca"], meta, [args.harmony_key], n_cells=adata.n_obs
    )
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)

    adata.write_h5ad(args.out, compression="gzip")
    U.log(f"RNA embed done: {adata.n_obs} cells, "
          f"{int(adata.var['highly_variable'].sum())} HVGs")


if __name__ == "__main__":
    main()
