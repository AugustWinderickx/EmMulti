#!/usr/bin/env python3
"""
atac_embed.py — filter + embed the ATAC modality.

Steps:
  * read the per-sample tile-matrix .h5ads produced by atac_import,
  * apply the per-sample TSSe threshold from atac_tsse_threshold,
  * scrublet + doublet removal,
  * combine into a SnapATAC2 AnnDataSet, select features, spectral embedding,
  * Harmony batch correction over ``sample`` (shape-safe),
  * build a gene-activity matrix (for marker cross-checks on the ATAC side),
  * emit an AnnData (obs + X_spectral / X_spectral_harmony in obsm) plus the
    gene-activity AnnData, both keyed by ``<study>_<sample>:<barcode>``.

Spectral-on-tiles is used (SnapATAC2 native) rather than LSI-on-peaks: it needs
no prior peak calling and is sufficient for broad WNN cell typing. Peak calling
can be added as a downstream module later.
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import snapatac2 as snap

import mwnn_utils as U


def sample_of(path):
    # efremova_CRC01.h5ad -> efremova_CRC01
    return Path(path).name.replace(".h5ad", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad_glob", required=True)
    ap.add_argument("--thresholds", required=True, help="tsse_thresholds.csv")
    ap.add_argument("--genome", default="hg38", choices=["hg38", "mm10"])
    ap.add_argument("--n_features", type=int, default=50000)
    ap.add_argument("--n_comps", type=int, default=30)
    ap.add_argument("--harmony_key", default="sample")
    ap.add_argument("--out_atac", required=True)
    ap.add_argument("--out_gene_activity", required=True)
    ap.add_argument("--combined_h5ads", required=True)
    args = ap.parse_args()

    genome = {"hg38": snap.genome.hg38, "mm10": snap.genome.mm10}[args.genome]
    thr = pd.read_csv(args.thresholds).set_index("sample")["applied_threshold"]

    files = sorted(glob.glob(args.h5ad_glob))
    if not files:
        raise SystemExit(f"no ATAC h5ads matched {args.h5ad_glob}")
    names = [sample_of(f) for f in files]

    adatas = [snap.read(f) for f in files]

    # per-sample TSSe filtering
    for name, a in zip(names, adatas):
        t = float(thr.get(name, thr.median()))
        snap.pp.filter_cells(a, min_tsse=t)
        U.log(f"{name}: TSSe>={t:.2f} -> {a.n_obs} cells")

    # feature selection + doublets per sample
    snap.pp.select_features(adatas, n_features=args.n_features)
    snap.pp.scrublet(adatas)
    snap.pp.filter_doublets(adatas)

    data = snap.AnnDataSet(
        adatas=list(zip(names, adatas)),
        filename=args.combined_h5ads,
        add_key="sample",
    )
    # obs_names already carry <study>_<sample>: prefix from import; guarantee unique
    data.obs_names = [f"{s}:{bc.split(':')[-1]}"
                      for s, bc in zip(data.obs["sample"], data.obs_names)]
    assert data.n_obs == np.unique(np.asarray(data.obs_names)).size

    snap.pp.select_features(data, n_features=args.n_features)
    snap.tl.spectral(data, n_comps=args.n_comps)

    # Harmony (shape-safe) on the spectral embedding
    meta = pd.DataFrame({args.harmony_key: np.asarray(data.obs[args.harmony_key])},
                        index=np.asarray(data.obs_names))
    Z = U.harmony_embedding(
        data.obsm["X_spectral"], meta, [args.harmony_key], n_cells=data.n_obs
    )
    data.obsm["X_spectral_harmony"] = Z

    # Build a plain AnnData carrying just what WNN needs (embeddings + obs)
    out = ad.AnnData(
        X=None,
        obs=pd.DataFrame(
            {"sample": np.asarray(data.obs["sample"])},
            index=np.asarray(data.obs_names),
        ),
        obsm={
            "X_spectral": np.asarray(data.obsm["X_spectral"]),
            "X_spectral_harmony": np.asarray(data.obsm["X_spectral_harmony"]),
        },
    )
    out.write_h5ad(args.out_atac)

    # Gene-activity matrix for ATAC-side marker cross-checks
    ga = snap.pp.make_gene_matrix(data, gene_anno=genome)
    ga.obs["sample"] = np.asarray(data.obs["sample"])
    import scanpy as sc
    sc.pp.normalize_total(ga)
    sc.pp.log1p(ga)
    ga.write_h5ad(args.out_gene_activity)

    U.log(f"ATAC embed done: {out.n_obs} cells, "
          f"spectral dims={out.obsm['X_spectral'].shape[1]}")
    data.close()
    for a in adatas:
        try:
            a.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
