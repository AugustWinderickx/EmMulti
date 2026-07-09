#!/usr/bin/env python3
"""
rna_import.py — read one sample's 10x .h5 into an AnnData with standardised
obs_names (<study>_<sample>:<barcode>) so it matches the ATAC side.
"""
import argparse

import numpy as np
import pandas as pd
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", required=True)
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--strip_barcode_suffix", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    a = sc.read_10x_h5(args.h5)
    a.var_names_make_unique()

    bc = pd.Series(np.asarray(a.obs_names, dtype=object))
    if args.strip_barcode_suffix:
        bc = bc.str.replace(r"-\d+$", "", regex=True)
    a.obs_names = np.array([f"{args.sample_id}:{b}" for b in bc])
    a.obs["sample"] = args.sample_id
    a.obs_names_make_unique()

    a.write_h5ad(args.out, compression="gzip")
    print(f"[rna_import] {args.sample_id}: {a.n_obs} cells x {a.n_vars} genes")


if __name__ == "__main__":
    main()
