#!/usr/bin/env python3
"""
rna_import.py — read one sample's RNA matrix into an AnnData with
standardised obs_names (<study>_<sample>:<barcode>) so it matches the ATAC
side.

Format-adaptive: most samplesheets point at a genuine CellRanger 10x .h5, but
external datasets sometimes ship an AnnData .h5ad (just not always named
.h5ad -- the Hickey/HuBMAP multiome release uses ``*_rna.h5``). The two have
incompatible internal layouts (10x's /matrix vs. AnnData's /X, /obs, /var, ...),
so the file is peeked with h5py first and read accordingly -- no extension
sniffing, no separate code path to remember to use for external samples.
"""
import argparse

import h5py
import numpy as np
import pandas as pd
import scanpy as sc

import mwnn_utils as U


def detect_format(path):
    with h5py.File(path, "r") as f:
        keys = set(f.keys())
    return "h5ad" if {"X", "obs", "var"}.issubset(keys) else "10x_h5"


def symbolize_if_needed(adata, gtf_path=None):
    """var_names -> gene symbol, for external AnnData inputs only.

    rna_qc_embed.py detects mt/ribo/hb genes and annotate.py scores marker
    sets by string-matching gene *symbols*. Some external .h5ad releases
    (e.g. Hickey/HuBMAP) carry versioned Ensembl IDs as var_names with the
    symbol in var['hugo_symbol'] instead; others carry bare/versioned Ensembl
    IDs with no symbol column at all -- against Ensembl IDs, all of those
    would silently match zero genes. No-op when var_names already look like
    symbols (the normal 10x_h5 case, and most external h5ads).

    Symbols are taken from var['hugo_symbol'] where present; anything still
    unresolved (no hugo_symbol column at all, or blank entries within it) is
    filled in from the GENCODE GTF at ``gtf_path``, keyed on the
    version-stripped Ensembl ID. Genes absent from both stay on their
    Ensembl ID rather than being dropped.
    """
    looks_ensembl = adata.var_names.str.match(r"^ENSG\d+(\.\d+)?$").mean() > 0.5
    if not looks_ensembl:
        return adata
    ens = adata.var_names.str.replace(r"\.\d+$", "", regex=True)
    if "hugo_symbol" in adata.var.columns:
        sym = adata.var["hugo_symbol"].astype(object)
    else:
        sym = pd.Series([None] * adata.n_vars, index=adata.var_names, dtype=object)
    missing = sym.isna() | (sym.astype(str).str.strip() == "")
    n_from_hugo = (~missing).sum()

    n_from_gtf = 0
    if missing.any() and gtf_path:
        gene_map = U.load_gtf_gene_map(gtf_path)
        filled = ens[missing].map(gene_map)
        sym.loc[missing] = filled.to_numpy()
        n_from_gtf = filled.notna().sum()
        missing = sym.isna() | (sym.astype(str).str.strip() == "")

    if missing.any() and not gtf_path:
        print(f"[rna_import] WARNING: {missing.sum()} Ensembl IDs unresolved "
              f"and no --gtf given -- falling back to Ensembl IDs for those genes")

    print(f"[rna_import] gene symbols: {n_from_hugo} from hugo_symbol, "
          f"{n_from_gtf} from GTF, {missing.sum()} unmapped (kept as Ensembl ID)")

    new_names = sym.where(~missing, ens)
    adata.var["ensembl_id"] = ens
    # .to_numpy() drops the Series name (would otherwise inherit "hugo_symbol"
    # and collide with the actual hugo_symbol column at write time).
    adata.var_names = pd.Index(new_names.astype(str).to_numpy())
    return adata


def read_rna(path, gtf_path=None):
    fmt = detect_format(path)
    a = sc.read_10x_h5(path) if fmt == "10x_h5" else symbolize_if_needed(sc.read_h5ad(path), gtf_path=gtf_path)
    a.var_names_make_unique()
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", required=True)
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--strip_barcode_suffix", action="store_true")
    ap.add_argument("--gtf", default=None,
                     help="GENCODE/Ensembl GTF used to fill in gene symbols "
                          "for inputs whose var_names are bare/versioned "
                          "Ensembl IDs (e.g. some external .h5ad releases) "
                          "-- covers genes missing from var['hugo_symbol'] "
                          "too, or files that don't carry that column at all.")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    a = read_rna(args.h5, gtf_path=args.gtf)
    U.assert_raw_counts(a, args.h5, check_marker=False)

    # rsplit on ":" so any pre-existing "<sample>:" prefix (external h5ads
    # that already went through their own import step) is dropped rather
    # than stacked; plain 10x barcodes have no ":" so this is a no-op there.
    bc = pd.Series(np.asarray(a.obs_names, dtype=object)).str.rsplit(":", n=1).str[-1]
    if args.strip_barcode_suffix:
        bc = bc.str.replace(r"-\d+$", "", regex=True)
    a.obs_names = np.array([f"{args.sample_id}:{b}" for b in bc])
    a.obs["sample"] = args.sample_id
    a.obs_names_make_unique()
    a.uns["emmulti_raw"] = True

    a.write_h5ad(args.out, compression="gzip")
    print(f"[rna_import] {args.sample_id}: {a.n_obs} cells x {a.n_vars} genes")


if __name__ == "__main__":
    main()
