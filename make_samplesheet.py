#!/usr/bin/env python3
"""
make_samplesheet.py — build a multiome-wnn-nf samplesheet from a
    <data_dir>/atac/<study>_<sample>_fragments.tsv.gz
    <data_dir>/rna/<study>_<sample>_rna.h5
layout, pairing samples and flagging anything unmatched or missing its index.

Usage:
    ./make_samplesheet.py --data_dir /…/multiome/efremova --out samplesheet.csv
    ./make_samplesheet.py --data_dir /…/multiome/efremova --study efremova \\
        --out samplesheet.csv
"""
import argparse
import sys
from pathlib import Path


ATAC_SUFFIX = "_fragments.tsv.gz"
RNA_SUFFIX = "_rna.h5"


def sample_from(name, suffix, study):
    stem = name[: -len(suffix)]              # efremova_CRC01
    prefix = f"{study}_"
    return stem[len(prefix):] if stem.startswith(prefix) else stem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True,
                    help="folder containing atac/ and rna/ subfolders")
    ap.add_argument("--study", default=None,
                    help="study name (default: data_dir folder name)")
    ap.add_argument("--out", default="samplesheet.csv")
    args = ap.parse_args()

    data_dir = Path(args.data_dir).resolve()
    study = args.study or data_dir.name
    atac_dir, rna_dir = data_dir / "atac", data_dir / "rna"
    for d in (atac_dir, rna_dir):
        if not d.is_dir():
            sys.exit(f"missing expected subfolder: {d}")

    atac = {sample_from(f.name, ATAC_SUFFIX, study): f
            for f in atac_dir.glob(f"*{ATAC_SUFFIX}")}
    rna = {sample_from(f.name, RNA_SUFFIX, study): f
           for f in rna_dir.glob(f"*{RNA_SUFFIX}")}

    paired = sorted(set(atac) & set(rna))
    atac_only = sorted(set(atac) - set(rna))
    rna_only = sorted(set(rna) - set(atac))

    # check tabix indexes
    missing_tbi = [s for s in paired
                   if not Path(str(atac[s]) + ".tbi").exists()]

    print(f"study            : {study}")
    print(f"paired samples   : {len(paired)}")
    if atac_only:
        print(f"ATAC-only (no RNA)  : {atac_only}")
    if rna_only:
        print(f"RNA-only (no ATAC)  : {rna_only}")
    if missing_tbi:
        print(f"WARNING missing .tbi: {missing_tbi}")

    if not paired:
        sys.exit("no paired samples found — check --study and file naming")

    with open(args.out, "w") as fh:
        fh.write("study,sample,atac_fragments,rna_h5\n")
        for s in paired:
            fh.write(f"{study},{s},{atac[s]},{rna[s]}\n")

    print(f"wrote {len(paired)} rows -> {args.out}")
    if missing_tbi:
        sys.exit("resolve missing .tbi files before running "
                 "(e.g. tabix -p bed <fragments>)")


if __name__ == "__main__":
    main()
