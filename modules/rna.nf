process RNA_IMPORT {
    tag { meta.id }
    label 'process_low'
    conda params.rna_env

    input:
    tuple val(meta), path(rna_h5)
    path gtf

    output:
    path "${meta.id}.rna.h5ad", emit: h5ad

    script:
    def strip = params.strip_barcode_suffix ? '--strip_barcode_suffix' : ''
    """
    rna_import.py \\
        --h5 ${rna_h5} \\
        --sample_id ${meta.id} \\
        ${strip} \\
        --gtf ${gtf} \\
        --out ${meta.id}.rna.h5ad
    """
}

process RNA_QC_EMBED {
    label 'process_high'
    conda params.rna_env
    publishDir "${params.outdir}/rna", mode: 'copy', pattern: '*.h5ad'
    publishDir "${params.outdir}/rna/qc", mode: 'copy', pattern: '*.png'

    input:
    path h5ads

    output:
    path "rna_embedded.h5ad", emit: rna
    path "*.png"

    script:
    def maxmt = params.rna_max_pct_mt ? "--max_pct_mt ${params.rna_max_pct_mt}" : ''
    """
    rna_qc_embed.py \\
        --h5ad_glob '*.rna.h5ad' \\
        --species ${params.species} \\
        --min_genes ${params.rna_min_genes} \\
        --min_cells ${params.rna_min_cells} \\
        ${maxmt} \\
        --qc_mode ${params.rna_qc_mode} \\
        --strict_n_mads ${params.rna_strict_n_mads} \\
        --strict_mt_n_mads ${params.rna_strict_mt_n_mads} \\
        --strict_max_pct_mt ${params.rna_strict_max_pct_mt} \\
        --n_top_genes ${params.rna_n_top_genes} \\
        --n_pcs ${params.rna_n_pcs} \\
        --outdir . \\
        --out rna_embedded.h5ad
    """
}
