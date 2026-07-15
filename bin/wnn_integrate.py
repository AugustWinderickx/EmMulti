#!/usr/bin/env python3
"""
wnn_integrate.py — join RNA + ATAC into a MuData and run WNN.

Because RNA and ATAC come from the same nuclei, cells are matched on
``<study>_<sample>:<barcode>``. Barcode-suffix mismatches (the classic ``-1``
pain) are normalised on both sides before intersecting.

Some 10x Multiome releases store fragments with the ATAC barcode as the
reverse complement of the barcode used to build the GEX/RNA library (a real,
observed quirk -- confirmed on Hickey/HuBMAP data, where direct matching gave
0% overlap for every sample). When direct matching is poor for a sample,
revcomp + the 10x ARC-v1 ATAC->GEX whitelist translation is tried instead.

Translation success (a barcode landing in the whitelist) does NOT by itself
mean the sample is genuinely cell-paired -- some samples in mixed cohorts are
separate, unpaired snRNA/snATAC captures on the same tissue block, not a true
single-reaction multiome. That shows up as post-translation overlap no better
than the *chance* rate for two random barcode sets of that size drawn from
the ~736k-barcode whitelist (matched, confirmed on Hickey B006-A-201: 31
observed vs. 27.8 expected by chance -- indistinguishable from noise -- vs.
sibling sample B006-A-001, same donor: 2397 observed vs. 35.5 expected, ~68x
enrichment, clearly real). So each sample's fate is decided per-sample by
fold-enrichment over chance, not by a hardcoded per-donor/per-study list --
a per-sample overlap report is written either way so silent losses are
visible and auditable.

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


def translate_all(obs_names, translation):
    """Per-cell revcomp + whitelist lookup; ``None`` where untranslatable."""
    out = np.empty(len(obs_names), dtype=object)
    for i, nm in enumerate(obs_names):
        pref, bc = nm.rsplit(":", 1)
        gex = translation.get(U.revcomp(bc))
        out[i] = f"{pref}:{gex}" if gex is not None else None
    return out


def resolve_pairing(rna, atac, translation, barcode_translation,
                     min_overlap_frac, min_fold):
    """Per-sample: decide direct / translated / unpaired, and return the
    (possibly barcode-translated) ATAC obs_names array plus an audit table.
    """
    atac_names = np.asarray(atac.obs_names, dtype=object)
    sample_arr = np.asarray(atac.obs["sample"])
    final_names = atac_names.copy()

    translated_all = None
    if translation is not None and barcode_translation != "skip":
        translated_all = translate_all(atac_names, translation)

    rep = []
    for s in sorted(set(rna.obs["sample"]) | set(atac.obs["sample"])):
        r = set(rna.obs_names[rna.obs["sample"] == s])
        smask = sample_arr == s
        a_direct = set(atac_names[smask])
        direct = len(r & a_direct)
        frac_direct = direct / max(len(r), 1)

        mode, matched, expected, fold = "direct", direct, None, None
        try_translate = translated_all is not None and (
            barcode_translation == "force" or frac_direct < min_overlap_frac
        )
        if try_translate:
            t_names = translated_all[smask]
            valid = t_names != None  # noqa: E711 (elementwise, not `is not None`)
            t_matched = len(r & set(t_names[valid]))
            n_translatable = int(valid.sum())
            expected = (n_translatable * len(r) / len(translation)
                        if translation else 0.0)
            fold = (t_matched / expected if expected > 0
                    else (float("inf") if t_matched > 0 else 0.0))
            if t_matched > direct and fold >= min_fold:
                mode, matched = "translated", t_matched
                sub = final_names[smask]
                sub[valid] = t_names[valid]
                final_names[smask] = sub
            elif direct == 0:
                mode, matched = "unpaired", 0

        rep.append({
            "sample": s, "rna": len(r), "atac": int(smask.sum()),
            "matched": matched, "mode": mode,
            "frac_rna": matched / max(len(r), 1),
            "frac_atac": matched / max(int(smask.sum()), 1),
            "expected_by_chance": None if expected is None else round(expected, 1),
            "fold_enrichment": None if fold is None else round(fold, 1),
        })
    return final_names, pd.DataFrame(rep)


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
                    help="warn if a sample matches below this fraction; also "
                         "the trigger, in 'auto' mode, for attempting "
                         "barcode translation on that sample")
    ap.add_argument("--barcode_translation", default="skip",
                    choices=["auto", "force", "skip"],
                    help="auto: try translation only where direct matching is "
                         "poor; force: always try; skip: never (legacy)")
    ap.add_argument("--atac_to_gex_translation", default=None,
                    help="gzipped TSV, 10x ARC-v1 ATAC->GEX barcode whitelist "
                         "translation (see assets/arc_whitelists/)")
    ap.add_argument("--barcode_pairing_min_fold", type=float, default=5.0,
                    help="min fold-enrichment of post-translation overlap "
                         "over the chance-level expectation to call a "
                         "sample genuinely cell-paired, vs. a separate/"
                         "unpaired snRNA+snATAC capture on the same tissue")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    rna = sc.read_h5ad(args.rna)
    atac = sc.read_h5ad(args.atac)

    rna.obs_names = normalise_names(rna.obs_names, args.strip_barcode_suffix)
    atac.obs_names = normalise_names(atac.obs_names, args.strip_barcode_suffix)

    translation = None
    if args.barcode_translation != "skip" and args.atac_to_gex_translation:
        translation = U.load_atac_gex_translation(args.atac_to_gex_translation)

    final_atac_names, rep = resolve_pairing(
        rna, atac, translation, args.barcode_translation,
        args.min_overlap_frac, args.barcode_pairing_min_fold,
    )
    atac.obs_names = pd.Index(final_atac_names)
    rep.to_csv(out / "wnn_overlap_report.csv", index=False)

    low = rep[(rep["frac_rna"] < args.min_overlap_frac) &
              (rep["rna"] > 0) & (rep["atac"] > 0)]
    for _, row in low.iterrows():
        if row["mode"] == "unpaired":
            U.log(f"WARNING {row['sample']}: unpaired -- 0 matched even after "
                  f"barcode translation ({row['fold_enrichment']}x chance "
                  f"expectation, need >={args.barcode_pairing_min_fold}x). "
                  f"Likely a separate snRNA/snATAC capture, not true "
                  f"multiome -- excluded from WNN.")
        else:
            U.log(f"WARNING low overlap for {row['sample']} ({row['mode']}): "
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
