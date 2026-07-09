#!/usr/bin/env python3
"""
atac_import.py — import one sample's ATAC fragments with SnapATAC2.

Per sample:
  * import_fragments (backed .h5ad)
  * compute per-cell TSS enrichment
  * add a 5 kb tile matrix (needed later for spectral)
  * export a small QC table (barcode, sample, n_fragments, tsse) that the
    threshold step aggregates across all samples.

Naming: obs_names are written as ``<study>_<sample>:<barcode>`` so that RNA and
ATAC line up cell-for-cell downstream.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import snapatac2 as snap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--sample_id", required=True, help="e.g. efremova_CRC01")
    ap.add_argument("--genome", default="hg38", choices=["hg38", "mm10"])
    ap.add_argument("--min_num_fragments", type=int, default=1000)
    ap.add_argument("--bin_size", type=int, default=5000)
    ap.add_argument("--strip_barcode_suffix", action="store_true",
                    help="strip trailing -<n> from barcodes for RNA/ATAC matching")
    ap.add_argument("--out_h5ad", required=True)
    ap.add_argument("--out_qc", required=True)
    args = ap.parse_args()

    genome = {"hg38": snap.genome.hg38, "mm10": snap.genome.mm10}[args.genome]

    adata = snap.pp.import_fragments(
        args.fragments,
        chrom_sizes=genome,
        min_num_fragments=args.min_num_fragments,
        file=args.out_h5ad,
        sorted_by_barcode=False,
    )

    snap.metrics.tsse(adata, genome)
    snap.pp.add_tile_matrix(adata, bin_size=args.bin_size)

    # standardise obs_names -> <study>_<sample>:<barcode>
    bc = np.asarray(adata.obs_names, dtype=object)
    if args.strip_barcode_suffix:
        bc = pd.Series(bc).str.replace(r"-\d+$", "", regex=True).to_numpy()
    adata.obs_names = np.array([f"{args.sample_id}:{b}" for b in bc])
    # snapatac2's backed AnnData (file=...) doesn't broadcast scalars like
    # pandas does; it needs an array matching adata.n_obs.
    adata.obs["sample"] = np.full(adata.n_obs, args.sample_id, dtype=object)

    # adata.obs is a Rust-backed PyDataFrameElem (polars), not pandas: no
    # .columns attribute, but __contains__ works directly.
    n_frag_col = "n_fragment" if "n_fragment" in adata.obs else None
    qc = pd.DataFrame({
        "cell_id": np.asarray(adata.obs_names),
        "sample": args.sample_id,
        "tsse": np.asarray(adata.obs["tsse"], dtype=float),
    })
    if n_frag_col:
        qc["n_fragments"] = np.asarray(adata.obs[n_frag_col], dtype=float)
    qc.to_parquet(args.out_qc, index=False)

    print(f"[atac_import] {args.sample_id}: {adata.n_obs} cells "
          f"(min_frag={args.min_num_fragments}), tsse median "
          f"{np.median(qc['tsse']):.2f}")
    adata.close()


if __name__ == "__main__":
    main()
