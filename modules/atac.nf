process ATAC_IMPORT {
    tag { meta.id }
    label 'process_medium'
    conda params.atac_env
    publishDir "${params.outdir}/atac/per_sample", mode: 'copy', pattern: '*.parquet'

    input:
    tuple val(meta), path(fragments), path(tbi)

    output:
    path "${meta.id}.h5ad",         emit: h5ad
    path "${meta.id}.qc.parquet",   emit: qc

    script:
    def strip = params.strip_barcode_suffix ? '--strip_barcode_suffix' : ''
    """
    atac_import.py \\
        --fragments ${fragments} \\
        --sample_id ${meta.id} \\
        --genome ${params.genome} \\
        --min_num_fragments ${params.min_num_fragments} \\
        --bin_size ${params.tile_size} \\
        ${strip} \\
        --out_h5ad ${meta.id}.h5ad \\
        --out_qc ${meta.id}.qc.parquet
    """
}

process ATAC_TSSE_THRESHOLD {
    label 'process_low'
    conda params.atac_env
    publishDir "${params.outdir}/atac/qc", mode: 'copy'

    input:
    path qc_files

    output:
    path "tsse_thresholds.csv", emit: thresholds
    path "*.png"
    path "tsse_threshold_meta.json"

    script:
    """
    atac_tsse_threshold.py \\
        --qc_glob '*.qc.parquet' \\
        --mode ${params.tsse_mode} \\
        --default_min_tsse ${params.default_min_tsse} \\
        --min_high_peak ${params.tsse_min_high_peak} \\
        --clamp_lo ${params.tsse_clamp_lo} \\
        --clamp_hi ${params.tsse_clamp_hi} \\
        --outdir .
    """
}

process ATAC_EMBED {
    label 'process_high'
    conda params.atac_env
    publishDir "${params.outdir}/atac", mode: 'copy', pattern: '*.h5ad'

    input:
    path h5ads
    path thresholds

    output:
    path "atac_embedded.h5ad",   emit: atac
    path "gene_activity.h5ad",   emit: gene_activity

    script:
    """
    atac_embed.py \\
        --h5ad_glob '*.h5ad' \\
        --thresholds ${thresholds} \\
        --genome ${params.genome} \\
        --n_features ${params.atac_n_features} \\
        --n_comps ${params.atac_n_comps} \\
        --out_atac atac_embedded.h5ad \\
        --out_gene_activity gene_activity.h5ad \\
        --combined_h5ads combined.h5ads
    """
}
