process WNN_INTEGRATE {
    label 'process_high'
    conda params.muon_env
    publishDir "${params.outdir}/wnn", mode: 'copy', pattern: '*.h5mu'
    publishDir "${params.outdir}/wnn/qc", mode: 'copy', pattern: '*.csv'

    input:
    path rna
    path atac

    output:
    path "mdata_wnn.h5mu", emit: mdata
    path "wnn_overlap_report.csv"

    script:
    def strip = params.strip_barcode_suffix ? '--strip_barcode_suffix' : ''
    """
    wnn_integrate.py \\
        --rna ${rna} \\
        --atac ${atac} \\
        --rna_rep ${params.rna_rep} \\
        --atac_rep ${params.atac_rep} \\
        --n_neighbors ${params.wnn_n_neighbors} \\
        --n_multineighbors ${params.wnn_n_multineighbors} \\
        ${strip} \\
        --outdir . \\
        --out mdata_wnn.h5mu
    """
}

process ANNOTATE {
    label 'process_medium'
    conda params.muon_env
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path mdata
    path markers
    path gene_activity

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
    label 'process_medium'
    conda params.muon_env
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path mdata
    path markers
    path gene_activity

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
