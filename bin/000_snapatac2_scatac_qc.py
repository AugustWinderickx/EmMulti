#!/usr/bin/env python3
import argparse
import snapatac2 as snap
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragments", required=True)
    parser.add_argument("--study_id", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--output_h5ad", required=True)
    parser.add_argument("--output_txt", required=True)
    args = parser.parse_args()

    adata = snap.pp.import_fragments(
        args.fragments, chrom_sizes=snap.genome.hg38,
        min_num_fragments=1000, file=args.output_h5ad, sorted_by_barcode=False
    )
    snap.metrics.tsse(adata, snap.genome.hg38)
    
    # Save TSSe array to a text file for global aggregation
    np.savetxt(args.output_txt, np.array(adata.obs["tsse"]))
    adata.close()

if __name__ == "__main__":
    main()