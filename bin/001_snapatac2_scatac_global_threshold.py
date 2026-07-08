#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs='+', required=True, help="List of TSSe text files")
    parser.add_argument("--study_id", required=True)
    parser.add_argument("--fallback_thresh", type=float, default=5.0)
    args = parser.parse_args()

    # Load and concatenate all TSSe scores
    all_tsse = np.concatenate([np.loadtxt(f) for f in args.inputs])
    
    # Calculate global dip
    data = all_tsse[(all_tsse >= 1) & (all_tsse <= 15)]
    global_thresh = args.fallback_thresh
    
    if len(data) >= 100:
        kde = gaussian_kde(data, bw_method=0.2)
        x_grid = np.linspace(1, 15, 1000)
        density_curve = kde(x_grid)
        peaks, _ = find_peaks(-density_curve, prominence=0.005)
        if len(peaks) > 0:
            global_thresh = round(x_grid[peaks[np.argmin(density_curve[peaks])]], 2)

    print(f"GLOBAL_THRESHOLD={global_thresh}")
    # Write just the number to a file for Nextflow to ingest
    with open("threshold.txt", "w") as f:
        f.write(str(global_thresh))

    # Generate and save the beautiful cohort-wide plot
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(all_tsse, bins=150, range=(0, 30), edgecolor="black", linewidth=0.2, color="#2c3e50", alpha=0.85)
    ax.axvline(global_thresh, color="#e74c3c", ls="--", lw=2, label=f"Global Dynamic Cutoff = {global_thresh}")
    ax.set_xlabel("TSSe Score")
    ax.set_ylabel("Cell Count")
    ax.set_title(f"Cohort TSSe Distribution & Global Threshold — Study: {args.study_id}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{args.study_id}_global_tsse_histogram.png", dpi=150)

if __name__ == "__main__":
    main()