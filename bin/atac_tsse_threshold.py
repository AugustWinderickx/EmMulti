#!/usr/bin/env python3
"""
atac_tsse_threshold.py — automatic TSS-enrichment thresholding.

Reads every per-sample QC parquet, then:
  * finds the antimode ("dip") of each sample's TSSe distribution,
  * finds a global dip on the pooled distribution,
  * writes a thresholds table (per-sample threshold, detected vs fallback),
  * draws adaptive QC figures whose size/rotation scale with #samples:
      - pooled TSSe histogram with the global dip marked
      - per-sample KDE ridge with each detected dip marked
      - per-sample TSSe boxplot with thresholds overlaid
      - (if n_fragments present) TSSe vs log10 fragments density knee plots
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mwnn_utils as U


def load_qc(qc_glob):
    files = sorted(glob.glob(qc_glob))
    if not files:
        raise SystemExit(f"no QC parquet files matched: {qc_glob}")
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    U.log(f"loaded {len(df):,} cells from {len(files)} samples")
    return df


def compute_thresholds(df, args):
    samples = sorted(df["sample"].unique())
    rows, per_sample_kde = [], {}
    for s in samples:
        v = df.loc[df["sample"] == s, "tsse"].to_numpy()
        r = U.find_tsse_threshold(
            v,
            search_range=(args.search_lo, args.search_hi),
            min_high_peak=args.min_high_peak,
            clamp=(args.clamp_lo, args.clamp_hi),
            default=args.default_min_tsse,
        )
        per_sample_kde[s] = r
        rows.append({
            "sample": s, "n_cells": r["n_cells"],
            "tsse_threshold": r["threshold"], "detected": r["detected"],
            "reason": r["reason"], "clamped": r["clamped"],
            "low_peak": r["low_peak"], "high_peak": r["high_peak"],
        })

    glob_r = U.find_tsse_threshold(
        df["tsse"].to_numpy(),
        search_range=(args.search_lo, args.search_hi),
        min_high_peak=args.min_high_peak,
        clamp=(args.clamp_lo, args.clamp_hi),
        default=args.default_min_tsse,
    )

    tbl = pd.DataFrame(rows)
    # Applied threshold depends on mode. In per_sample mode, fall back to the
    # global dip for any sample where auto-detection failed. In manual mode,
    # the user-supplied value overrides everywhere (auto-detection above is
    # still run and kept in the table for reference/diagnostics only).
    if args.mode == "manual":
        if args.manual_threshold is None:
            raise SystemExit("--manual_threshold is required when --mode manual")
        tbl["applied_threshold"] = args.manual_threshold
    elif args.mode == "global":
        tbl["applied_threshold"] = glob_r["threshold"]
    else:
        tbl["applied_threshold"] = np.where(
            tbl["detected"], tbl["tsse_threshold"], glob_r["threshold"]
        )
    return tbl, per_sample_kde, glob_r


def plot_pooled_hist(df, threshold, label, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["tsse"], bins=120, range=(0, 30),
            color="#7089b8", edgecolor="white", linewidth=0.2)
    ax.axvline(threshold, color="#d1495b", ls="--", lw=2, label=label)
    ax.set(xlabel="TSS enrichment", ylabel="cells",
           title="Pooled TSSe distribution")
    ax.legend(frameon=False)
    U.save_fig(fig, path)


def plot_sample_ridge(kde, path):
    samples = list(kde.keys())
    n = len(samples)
    fig, ax = plt.subplots(figsize=(9, max(3, 0.42 * n + 1)))
    colors = U.categorical_palette(n)
    offset_step = 1.0
    for i, s in enumerate(samples):
        r = kde[s]
        if r["grid_y"] is None:
            continue
        x = np.asarray(r["grid_x"]); y = np.asarray(r["grid_y"])
        y = y / (y.max() + 1e-12)  # normalise height
        base = i * offset_step
        ax.fill_between(x, base, base + y, color=colors[i], alpha=0.65, lw=0)
        ax.plot(x, base + y, color="black", lw=0.4)
        ax.plot([r["threshold"]], [base], marker="v", color="#d1495b", ms=7)
    ax.set_yticks([i * offset_step for i in range(n)])
    ax.set_yticklabels(samples, fontsize=max(5, min(9, int(9 - n / 40))))
    ax.set(xlabel="TSS enrichment",
           title="Per-sample TSSe (▼ = detected dip)")
    ax.set_xlim(0, None)
    U.save_fig(fig, path)


def plot_sample_boxplot(df, tbl, path):
    samples = list(tbl["sample"])
    n = len(samples)
    data = [df.loc[df["sample"] == s, "tsse"].to_numpy() for s in samples]
    fig, ax = plt.subplots(figsize=(U.sample_figwidth(n), 5))
    ax.boxplot(data, showfliers=False, widths=0.6)
    thr = tbl.set_index("sample")["applied_threshold"]
    ax.scatter(range(1, n + 1), thr.loc[samples].to_numpy(),
               color="#d1495b", zorder=5, s=22, label="applied threshold")
    U.style_sample_xticks(ax, samples)
    ax.set(ylabel="TSS enrichment", title="TSSe by sample")
    ax.legend(frameon=False)
    U.save_fig(fig, path)


def plot_knee(df, tbl, path):
    if "n_fragments" not in df.columns:
        return
    samples = list(tbl["sample"])
    nr, nc = U.grid_dims(len(samples), max_cols=4)
    fig, axes = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.0 * nr), squeeze=False)
    thr = tbl.set_index("sample")["applied_threshold"]
    for ax, s in zip(axes.ravel(), samples):
        d = df[df["sample"] == s]
        ax.scatter(np.log10(d["n_fragments"] + 1), d["tsse"], s=1.5,
                   alpha=0.25, color="#40587c", rasterized=True)
        ax.axhline(thr[s], color="#d1495b", ls="--", lw=1)
        ax.set_title(s, fontsize=8)
        ax.set_xlabel("log10 fragments", fontsize=7)
        ax.set_ylabel("TSSe", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[len(samples):]:
        ax.axis("off")
    fig.suptitle("TSSe vs fragment count (dashed = applied threshold)")
    U.save_fig(fig, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc_glob", required=True)
    ap.add_argument("--mode", default="per_sample",
                    choices=["per_sample", "global", "manual"])
    ap.add_argument("--manual_threshold", type=float, default=None,
                    help="fixed TSSe threshold applied to every sample "
                         "when --mode manual")
    ap.add_argument("--default_min_tsse", type=float, default=5.0)
    ap.add_argument("--min_high_peak", type=float, default=4.0)
    ap.add_argument("--search_lo", type=float, default=0.0)
    ap.add_argument("--search_hi", type=float, default=20.0)
    ap.add_argument("--clamp_lo", type=float, default=2.0)
    ap.add_argument("--clamp_hi", type=float, default=10.0)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    df = load_qc(args.qc_glob)
    tbl, kde, glob_r = compute_thresholds(df, args)

    tbl.to_csv(out / "tsse_thresholds.csv", index=False)
    U.write_json({"mode": args.mode, "global": glob_r,
                  "params": vars(args)}, out / "tsse_threshold_meta.json")

    if args.mode == "manual":
        pooled_thr, pooled_label = args.manual_threshold, \
            f"manual threshold = {args.manual_threshold:.2f}"
    else:
        pooled_thr = glob_r["threshold"]
        pooled_label = (f"global dip = {pooled_thr:.2f}"
                         + ("" if glob_r["detected"] else " (fallback)"))
    plot_pooled_hist(df, pooled_thr, pooled_label, out / "tsse_pooled_hist.png")
    plot_sample_ridge(kde, out / "tsse_sample_ridge.png")
    plot_sample_boxplot(df, tbl, out / "tsse_sample_boxplot.png")
    plot_knee(df, tbl, out / "tsse_knee.png")

    n_det = int(tbl["detected"].sum())
    U.log(f"thresholds: {n_det}/{len(tbl)} auto-detected; "
          f"global dip = {glob_r['threshold']:.2f} "
          f"({'detected' if glob_r['detected'] else glob_r['reason']})")


if __name__ == "__main__":
    main()
