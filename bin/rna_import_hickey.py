#!/usr/bin/env python3
"""
rna_import_hickey.py — import HuBMAP/Hickey healthy-gut scRNA *_raw_expr.h5ad
files into the same <sample_id>.rna.h5ad convention rna_import.py produces,
so RNA_QC_EMBED's ``*.rna.h5ad`` glob picks them up unchanged.

Each source file has:
  * var_names = versioned Ensembl IDs (e.g. ENSG00000000003.15), with the
    HGNC symbol in var['hugo_symbol'] -- missing for ~42% of genes
    (non-coding / unannotated), so a straight rename would produce NaN names,
  * bare-barcode obs_names (no ``-1`` suffix) and no obs columns at all.

rna_qc_embed.py detects mt/ribo/hb genes and annotate.py scores marker sets
by string-matching gene *symbols* (assets/markers.yaml), so var_names here
are remapped to hugo_symbol -- falling back to the (version-stripped)
Ensembl ID for genes with no symbol, so no gene is silently dropped.
"""
import argparse
import glob
from pathlib import Path

import pandas as pd
import scanpy as sc

import mwnn_utils as U


def symbolize(adata, symbol_map=None):
    """var_names -> hugo_symbol, falling back to the Ensembl ID when absent.

    A minority of Hickey files carry no ``var['hugo_symbol']`` column at all
    (7 of 78, checked directly). All files share the exact same var_names /
    gene order (same GENCODE v35 reference), so ``symbol_map`` -- built once
    from a file that does have the column -- covers those too, instead of
    silently falling back to bare Ensembl IDs for those samples only.
    """
    ens = adata.var_names.str.replace(r"\.\d+$", "", regex=True)
    if "hugo_symbol" in adata.var.columns:
        sym = adata.var["hugo_symbol"].astype(object)
    elif symbol_map is not None:
        sym = pd.Series(symbol_map).reindex(adata.var_names).astype(object)
    else:
        sym = pd.Series([None] * adata.n_vars, index=adata.var_names, dtype=object)
    missing = sym.isna() | (sym.astype(str).str.strip() == "")
    new_names = sym.where(~missing, ens)
    adata.var["hugo_symbol"] = sym.to_numpy()
    adata.var["ensembl_id"] = ens
    # .to_numpy() drops the Series name -- pd.Index(name=None) alone doesn't,
    # since None is also its unpassed default and it silently inherits the
    # Series' name ("hugo_symbol"), which then collides with the actual
    # hugo_symbol column at write time.
    adata.var_names = pd.Index(new_names.astype(str).to_numpy())
    adata.var_names_make_unique()
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", required=True,
                    help="e.g. '.../hickey/*/*_raw_expr.h5ad' (covers both "
                         "small_intestine/ and large_intestine/)")
    ap.add_argument("--study", default="hickey")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(args.glob))
    if not files:
        raise SystemExit(f"no files matched {args.glob}")

    # Build one Ensembl->symbol map from whichever file has var['hugo_symbol'],
    # to cover the files that don't (all files share the same reference, so
    # this is exact, not a guess).
    symbol_map = None
    for f in files:
        probe = sc.read_h5ad(f, backed="r")
        if "hugo_symbol" in probe.var.columns:
            symbol_map = probe.var["hugo_symbol"].astype(object)
            break
    if symbol_map is None:
        print("[rna_import_hickey] WARNING: no file has var['hugo_symbol']; "
              "falling back to Ensembl IDs everywhere")

    for f in files:
        f = Path(f)
        stem = f.name.removesuffix("_raw_expr.h5ad")
        tissue = f.parent.name  # small_intestine | large_intestine
        sample_id = f"{args.study}_{tissue}_{stem}"

        a = sc.read_h5ad(f)
        U.assert_raw_counts(a, str(f), check_marker=False)
        a = symbolize(a, symbol_map=symbol_map)
        a.obs["sample"] = sample_id
        a.obs["tissue"] = tissue
        a.obs_names = [f"{sample_id}:{bc}" for bc in a.obs_names]
        a.obs_names_make_unique()
        a.uns["emmulti_raw"] = True

        out_path = out / f"{sample_id}.rna.h5ad"
        a.write_h5ad(out_path, compression="gzip")
        print(f"[rna_import_hickey] {sample_id}: "
              f"{a.n_obs} cells x {a.n_vars} genes -> {out_path}")


if __name__ == "__main__":
    main()
