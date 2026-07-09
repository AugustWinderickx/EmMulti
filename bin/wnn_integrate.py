#!/usr/bin/env python3
"""
wnn_integrate.py — join RNA + ATAC into a MuData and run WNN.

Because RNA and ATAC come from the same nuclei, cells are matched on
``<study>_<sample>:<barcode>``. Barcode-suffix mismatches (the classic ``-1``
pain) are normalised on both sides before intersecting; a per-sample overlap
report is written so silent losses are visible.

Then:
  * per-modality neighbours (RNA: X_pca_harmony, ATAC: X_spectral_harmony),
  * mu.pp.neighbors (weighted nearest neighbours),
  * mu.tl.umap.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu

import mwnn_utils as U


def normalise_names(names, strip_suffix):
    s = pd.Series(np.asarray(names, dtype=object))
    if strip_suffix:
        # only touch the barcode part after the last ':'
        pref = s.str.rsplit(":", n=1).str[0]
        bc = s.str.rsplit(":", n=1).str[-1].str.replace(r"-\d+$", "", regex=True)
        s = pref + ":" + bc
    return s.to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", required=True)
    ap.add_argument("--atac", required=True)
    ap.add_argument("--rna_rep", default="X_pca_harmony")
    ap.add_argument("--atac_rep", default="X_spectral_harmony")
    ap.add_argument("--n_neighbors", type=int, default=50)
    ap.add_argument("--n_multineighbors", type=int, default=50)
    ap.add_argument("--strip_barcode_suffix", action="store_true")
    ap.add_argument("--min_overlap_frac", type=float, default=0.5,
                    help="warn if a sample matches below this fraction")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    rna = sc.read_h5ad(args.rna)
    atac = sc.read_h5ad(args.atac)

    rna.obs_names = normalise_names(rna.obs_names, args.strip_barcode_suffix)
    atac.obs_names = normalise_names(atac.obs_names, args.strip_barcode_suffix)

    # per-sample overlap report
    rep = []
    for s in sorted(set(rna.obs["sample"]) | set(atac.obs["sample"])):
        r = set(rna.obs_names[rna.obs["sample"] == s])
        a = set(atac.obs_names[atac.obs["sample"] == s])
        inter = len(r & a)
        rep.append({"sample": s, "rna": len(r), "atac": len(a),
                    "matched": inter,
                    "frac_rna": inter / max(len(r), 1),
                    "frac_atac": inter / max(len(a), 1)})
    rep = pd.DataFrame(rep)
    rep.to_csv(out / "wnn_overlap_report.csv", index=False)
    low = rep[(rep["frac_rna"] < args.min_overlap_frac) &
              (rep["rna"] > 0) & (rep["atac"] > 0)]
    for _, row in low.iterrows():
        U.log(f"WARNING low overlap for {row['sample']}: "
              f"{row['matched']} matched ({row['frac_rna']:.0%} of RNA)")

    common = rna.obs_names.intersection(atac.obs_names)
    if len(common) == 0:
        raise SystemExit("no cells matched between RNA and ATAC — check naming "
                         "and --strip_barcode_suffix")
    rna = rna[common].copy()
    atac = atac[common].copy()
    U.log(f"matched {len(common):,} multiome cells")

    mdata = mu.MuData({"rna": rna, "atac": atac})
    mdata.obs["sample"] = rna.obs.loc[mdata.obs_names, "sample"].astype(str).values

    sc.pp.neighbors(mdata["rna"], use_rep=args.rna_rep,
                    n_neighbors=args.n_neighbors, metric="cosine",
                    key_added="rna_neighbors")
    sc.pp.neighbors(mdata["atac"], use_rep=args.atac_rep,
                    n_neighbors=args.n_neighbors, metric="cosine",
                    key_added="atac_neighbors")
    mu.pp.neighbors(
        mdata,
        neighbor_keys={"rna": "rna_neighbors", "atac": "atac_neighbors"},
        n_neighbors=args.n_neighbors,
        n_multineighbors=args.n_multineighbors,
        metric="cosine",
    )
    mu.tl.umap(mdata, min_dist=0.1, spread=1.0,
               init_pos="spectral", random_state=0)

    mdata.write_h5mu(args.out)
    U.log(f"WNN done: {mdata.n_obs} cells")


if __name__ == "__main__":
    main()
