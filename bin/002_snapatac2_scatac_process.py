#!/usr/bin/env python3
import argparse
import snapatac2 as snap

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_h5ad", required=True)
    parser.add_argument("--global_thresh", type=float, required=True)
    parser.add_argument("--study_id", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--output_h5ad", required=True)
    args = parser.parse_args()

    # Reopen backed file
    adata = snap.read(args.input_h5ad, backed="r+")
    
    # Filter with global threshold
    snap.pp.filter_cells(adata, min_tsse=args.global_thresh)
    
    # Complete downstream processing
    snap.pp.add_tile_matrix(adata, bin_size=5000)
    snap.pp.scrublet(adata)
    snap.pp.filter_doublets(adata)
    
    # Standardize Metadata Columns
    adata.obs['study'] = args.study_id
    adata.obs['sample'] = args.sample_id
    adata.obs['raw_barcode'] = adata.obs_names
    adata.obs_names = (adata.obs["study"].astype(str) + ":" + 
                       adata.obs["sample"].astype(str) + ":" + 
                       adata.obs["raw_barcode"].astype(str))
    adata.close()

if __name__ == "__main__":
    main()