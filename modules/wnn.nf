process WNN_INTEGRATE {
    tag { variant }
    label 'process_high'
    conda params.muon_env
    publishDir "${params.outdir}/wnn/${variant}", mode: 'copy', pattern: '*.h5mu'
    publishDir "${params.outdir}/wnn/${variant}/qc", mode: 'copy', pattern: '*.csv'

    input:
    path rna
    path atac
    val variant
    val rna_rep
    val atac_rep
    path atac_to_gex_translation

    output:
    path "mdata_wnn.h5mu", emit: mdata
    path "wnn_overlap_report.csv"

    script:
    def strip = params.strip_barcode_suffix ? '--strip_barcode_suffix' : ''
    """
    wnn_integrate.py \\
        --rna ${rna} \\
        --atac ${atac} \\
        --rna_rep ${rna_rep} \\
        --atac_rep ${atac_rep} \\
        --n_neighbors ${params.wnn_n_neighbors} \\
        --n_multineighbors ${params.wnn_n_multineighbors} \\
        ${strip} \\
        --barcode_translation ${params.barcode_translation} \\
        --atac_to_gex_translation ${atac_to_gex_translation} \\
        --barcode_pairing_min_fold ${params.barcode_pairing_min_fold} \\
        --outdir . \\
        --out mdata_wnn.h5mu
    """
}

process ANNOTATE {
    tag { variant }
    label 'process_medium'
    conda params.muon_env
    publishDir "${params.outdir}/annotation/${variant}", mode: 'copy'

    input:
    path mdata
    path markers
    path gene_activity
    val variant

    output:
    path "mdata_annotated.h5mu", emit: mdata
    path "*.csv"
    path "cell_annotations.tsv"

    script:
    """
    annotate.py \\
        --mdata ${mdata} \\
        --markers ${markers} \\
        --gene_activity ${gene_activity} \\
        --resolutions ${params.leiden_resolutions} \\
        --primary_resolution ${params.primary_resolution} \\
        --outdir . \\
        --out mdata_annotated.h5mu
    """
}

process PLOTS {
    tag { variant }
    label 'process_medium'
    conda params.muon_env
    publishDir "${params.outdir}/plots/${variant}", mode: 'copy'

    input:
    path mdata
    path markers
    path gene_activity
    val variant

    output:
    path "*.png"

    script:
    """
    plots.py \\
        --mdata ${mdata} \\
        --markers ${markers} \\
        --gene_activity ${gene_activity} \\
        --primary_resolution ${params.primary_resolution} \\
        --outdir .
    """
}
