"""
mwnn_utils.py — shared helpers for the multiome-wnn-nf pipeline.

Grouped into:
  * IO / logging helpers
  * harmony safe-transpose (fixes the Z_corr shape gotcha)
  * adaptive plotting (figsize / grids / palettes / point sizes / legends
    that scale with the number of samples so text never overlaps)
  * TSSe auto-threshold via KDE antimode ("the dip") detection

Everything here is import-safe: no side effects at import time except a
non-interactive matplotlib backend.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")  # headless on the cluster
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex


# --------------------------------------------------------------------------- #
# IO / logging
# --------------------------------------------------------------------------- #
def log(*msg):
    print("[mwnn]", *msg, file=sys.stderr, flush=True)


def write_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default))


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# --------------------------------------------------------------------------- #
# raw-counts guard
# --------------------------------------------------------------------------- #
def assert_raw_counts(adata, source_name, check_marker=True):
    """Guard against double-normalizing input that isn't raw counts.

    normalize_total/log1p run unconditionally downstream (rna_qc_embed.py),
    with no way to tell after the fact whether they're being applied to raw
    counts or to already-normalized/scaled data someone pointed the
    samplesheet at by mistake (external .h5ad releases vary in what state
    they ship in -- see the Hickey secondary_analysis.h5ad case, whose X was
    already z-scored down to -3.07).

    Negative values are an unambiguous tell (raw counts are never negative)
    and are a hard failure. Missing the emmulti_raw marker -- stamped by
    rna_import.py / rna_import_hickey.py on the way in -- is only a warning,
    since legitimately-raw data can arrive without having gone through them.
    ``check_marker=False`` at import time, before the marker has been set.
    """
    xmin = adata.X.min()
    if xmin < 0:
        raise SystemExit(
            f"{source_name}: X has negative values (min={xmin:.3g}) -- looks "
            f"already scaled/z-scored, not raw counts. Refusing to "
            f"normalize_total/log1p on top of it."
        )
    if check_marker and not adata.uns.get("emmulti_raw", False):
        log(f"WARNING {source_name}: no 'emmulti_raw' marker in .uns -- this "
            f"file didn't come through rna_import.py/rna_import_hickey.py, "
            f"so raw-counts status can't be confirmed. Proceeding.")


# --------------------------------------------------------------------------- #
# 10x Multiome ATAC<->GEX barcode translation
# --------------------------------------------------------------------------- #
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def load_atac_gex_translation(path):
    """Load the 10x ARC-v1 ATAC->GEX barcode translation table.

    ``path`` is a gzipped two-column TSV (ATAC barcode, GEX barcode), built
    by pairing the two whitelists positionally -- see assets/arc_whitelists/.
    Fragments files store the ATAC barcode reverse-complemented relative to
    the whitelist (a real, observed quirk, not universal 10x behaviour), so
    callers must revcomp before looking up here.
    """
    import gzip
    table = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            atac_bc, gex_bc = line.rstrip("\n").split("\t")
            table[atac_bc] = gex_bc
    return table


# --------------------------------------------------------------------------- #
# harmony
# --------------------------------------------------------------------------- #
def harmony_embedding(X, meta_df, vars_use, n_cells=None, max_iter=20, seed=0, **kw):
    """
    Run harmonypy and ALWAYS return a (n_cells, n_dims) array.

    harmonypy's ``run_harmony(...).Z_corr`` is shaped (n_dims, n_cells), so
    assigning it straight into ``obsm`` (which expects rows == cells) either
    errors or silently corrupts the embedding. This wrapper transposes as
    needed and validates against the known cell count.
    """
    import harmonypy as hm

    X = np.asarray(X, dtype=np.float64)
    if n_cells is None:
        n_cells = X.shape[0]
    ho = hm.run_harmony(X, meta_df, vars_use, max_iter_harmony=max_iter,
                        random_state=seed, **kw)
    Z = np.asarray(ho.Z_corr)
    Z = _orient_to_cells(Z, n_cells)
    return Z


def _orient_to_cells(Z, n_cells):
    if Z.ndim != 2:
        raise ValueError(f"expected 2D harmony output, got shape {Z.shape}")
    if Z.shape[0] == n_cells:
        return np.ascontiguousarray(Z)
    if Z.shape[1] == n_cells:
        return np.ascontiguousarray(Z.T)
    raise ValueError(
        f"harmony output shape {Z.shape} matches n_cells={n_cells} on neither axis"
    )


# --------------------------------------------------------------------------- #
# Adaptive plotting
# --------------------------------------------------------------------------- #
def grid_dims(n, ncols=None, max_cols=4):
    """Rows/cols for n panels; keeps panels squarish and readable."""
    if n <= 0:
        return 1, 1
    if ncols is None:
        ncols = min(max_cols, n)
    nrows = int(np.ceil(n / ncols))
    return nrows, ncols


def categorical_palette(n):
    """
    Return n visually distinct hex colors. Uses colorcet's glasbey palette
    when available (best for many categories), otherwise stitches together
    matplotlib qualitative maps, finally falling back to evenly-spaced HSV.
    """
    n = int(n)
    if n <= 0:
        return []
    try:
        import colorcet as cc
        pal = cc.glasbey_bw_minc_20 if n > 20 else cc.glasbey_category10
        if len(pal) >= n:
            return [to_hex(c) for c in pal[:n]]
    except Exception:
        pass

    base = []
    for cmap_name in ("tab20", "tab20b", "tab20c", "Set3", "Dark2"):
        cmap = plt.get_cmap(cmap_name)
        base.extend(cmap.colors)
        if len(base) >= n:
            return [to_hex(c) for c in base[:n]]

    # last resort: evenly spaced hues
    import colorsys
    return [to_hex(colorsys.hsv_to_rgb(i / n, 0.65, 0.9)) for i in range(n)]


def sample_figwidth(n_samples, per_sample=0.42, base=4.0, lo=6.0, hi=42.0):
    """Width (inches) for a per-sample axis (boxplots, bar charts)."""
    return float(np.clip(base + per_sample * n_samples, lo, hi))


def style_sample_xticks(ax, labels, max_before_rotate=8):
    """Rotate + shrink x tick labels based on count and length so they never
    collide. Long/many labels -> 90 deg + smaller font."""
    n = len(labels)
    longest = max((len(str(x)) for x in labels), default=0)
    if n <= max_before_rotate and longest <= 10:
        rot, fs, ha = 0, 10, "center"
    elif n <= 24:
        rot, fs, ha = 45, 9, "right"
    else:
        rot, fs, ha = 90, max(5, int(8 - n / 40)), "center"
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=rot, fontsize=fs, ha=ha)


def legend_outside(ax, n_items, title=None):
    """Put a legend to the right; wrap into multiple columns when it's long."""
    ncol = 1 if n_items <= 22 else 2 if n_items <= 44 else 3
    fs = 9 if n_items <= 22 else 7 if n_items <= 60 else 6
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              ncol=ncol, fontsize=fs, frameon=False, title=title,
              markerscale=1.5, handletextpad=0.3, columnspacing=0.8)


def umap_point_size(n_cells):
    """Point size that keeps large embeddings from turning into a solid blob."""
    if n_cells <= 2_000:
        return 20.0
    if n_cells <= 20_000:
        return 8.0
    if n_cells <= 100_000:
        return 3.0
    if n_cells <= 400_000:
        return 1.2
    return 0.5


def save_fig(fig, path, dpi=200):
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log("wrote", path)


# --------------------------------------------------------------------------- #
# TSSe auto-threshold ("the dip")
# --------------------------------------------------------------------------- #
def find_tsse_threshold(
    tsse,
    search_range=(0.0, 20.0),
    grid=512,
    min_high_peak=4.0,
    clamp=(2.0, 10.0),
    default=5.0,
    min_cells=50,
    bw_method=None,
):
    """
    Locate the antimode ("dip") between the low-quality and high-quality modes
    of a per-cell TSS-enrichment distribution.

    Strategy
    --------
    1. Restrict to a sensible window (default 0-20) and KDE-smooth the values.
    2. Find KDE peaks (modes) and valleys (antimodes).
    3. Pick the dominant *low* peak (< ``min_high_peak``) and dominant *high*
       peak (>= ``min_high_peak``); the deepest valley between them is the
       threshold.
    4. Clamp into a safe range and fall back to ``default`` if the distribution
       is not clearly bimodal, always reporting which happened.

    Returns a dict with the threshold, whether it was auto-detected, the KDE
    grid + density (for plotting), and the peak/valley locations.
    """
    from scipy.stats import gaussian_kde
    from scipy.signal import find_peaks

    tsse = np.asarray(tsse, dtype=float)
    tsse = tsse[np.isfinite(tsse)]
    lo, hi = search_range
    x = np.linspace(lo, hi, grid)

    result = {
        "threshold": float(np.clip(default, *clamp)),
        "detected": False,
        "reason": "",
        "n_cells": int(tsse.size),
        "grid_x": x.tolist(),
        "grid_y": None,
        "low_peak": None,
        "high_peak": None,
        "valley": None,
        "clamped": False,
    }

    vals = tsse[(tsse >= lo) & (tsse <= hi)]
    if vals.size < min_cells:
        result["reason"] = "too_few_cells"
        return result

    # subsample for KDE speed on huge samples (deterministic)
    if vals.size > 200_000:
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, 200_000, replace=False)

    try:
        kde = gaussian_kde(vals, bw_method=bw_method)
        y = kde(x)
    except Exception as e:  # singular data etc.
        result["reason"] = f"kde_failed:{type(e).__name__}"
        return result
    result["grid_y"] = y.tolist()

    peak_idx, _ = find_peaks(y)
    if peak_idx.size < 2:
        result["reason"] = "not_bimodal"
        return result

    px, py = x[peak_idx], y[peak_idx]
    low_sel = px < min_high_peak
    high_sel = px >= min_high_peak
    if not low_sel.any() or not high_sel.any():
        result["reason"] = "no_low_and_high_mode"
        return result

    lp = peak_idx[low_sel][np.argmax(py[low_sel])]
    hp = peak_idx[high_sel][np.argmax(py[high_sel])]
    if hp <= lp:
        result["reason"] = "modes_out_of_order"
        return result

    seg = slice(lp, hp + 1)
    v = lp + int(np.argmin(y[seg]))
    thr = float(x[v])
    clamped = not (clamp[0] <= thr <= clamp[1])
    thr_clamped = float(np.clip(thr, *clamp))

    result.update(
        threshold=thr_clamped,
        detected=True,
        reason="ok",
        low_peak=float(x[lp]),
        high_peak=float(x[hp]),
        valley=thr,
        clamped=clamped,
    )
    return result
