nextflow.enable.dsl = 2

params.atac_dir   = "/lustre1/project/stg_00064/awinderickx/data/crc_reverse/metastasis/multiome/efremova/atac"
params.outdir     = "results"

process ATAC_STAGE1_QC {
    tag "${study}_${sample}"
    
    input:
    tuple val(study), val(sample), path(fragments)

    output:
    tuple val(study), val(sample), path("${study}_${sample}_raw.h5ad"), emit: h5ad
    path "${study}_${sample}_tsse.txt", emit: tsse_txt

    script:
    """
    atac_stage1_qc.py \\
        --fragments ${fragments} \\
        --study_id ${study} \\
        --sample_id ${sample} \\
        --output_h5ad ${study}_${sample}_raw.h5ad \\
        --output_txt ${study}_${sample}_tsse.txt
    """
}

process ATAC_STAGE2_GLOBAL_THRESHOLD {
    tag "${study}_global_qc"
    publishDir "${params.outdir}/${study}/atac/qc_plots", mode: 'copy', pattern: "*.png"
    
    input:
    tuple val(study), path(tsse_files)

    output:
    path "threshold.txt", emit: thresh_val
    path "${study}_global_tsse_histogram.png"

    script:
    """
    atac_stage2_global_threshold.py \\
        --inputs ${tsse_files} \\
        --study_id ${study}
    """
}

process ATAC_STAGE3_PROCESS {
    tag "${study}_${sample}"
    
    input:
    tuple val(study), val(sample), path(raw_h5ad)
    val global_thresh

    output:
    tuple val(study), path("${study}_${sample}_processed.h5ad")

    script:
    """
    atac_stage3_process.py \\
        --input_h5ad ${raw_h5ad} \\
        --global_thresh ${global_thresh} \\
        --study_id ${study} \\
        --sample_id ${sample} \\
        --output_h5ad ${study}_${sample}_processed.h5ad
    """
}

process MERGE_AND_HARMONY_ATAC {
    tag "${study}_atlas"
    publishDir "${params.outdir}/${study}/atac", mode: 'copy'

    input:
    tuple val(study), path(h5ad_files)

    output:
    path "${study}_combined_atac.h5ad"

    script:
    """
    merge_atac_and_embed.py \\
        --inputs ${h5ad_files} \\
        --output ${study}_combined_atac.h5ad
    """
}

workflow {
    atac_ch = Channel
        .fromPath("${params.atac_dir}/*_*_fragments.tsv.gz")
        .map { file ->
            def parts = file.name.tokenize('_')
            tuple(parts[0], parts[1], file) // study, sample, file
        }

    // 1. Calculate raw metric files in parallel across Slurm nodes
    stage1_out = ATAC_STAGE1_QC(atac_ch)

    // 2. Group all TSSe text files by study and calculate the global cohort threshold
    study_grouped_txt = stage1_out.tsse_txt
        .map { file -> 
            def study = file.name.tokenize('_')[0]
            tuple(study, file) 
        }.groupTuple()
        
    global_qc = ATAC_STAGE2_GLOBAL_THRESHOLD(study_grouped_txt)

    // Ingest the computed cutoff string as a dynamic pipeline variable
    threshold_value_ch = global_qc.thresh_val.map { it.text.trim() }

    // 3. Re-combine original raw h5ads with the global threshold for parallel filtering & feature processing
    // We mix the sample channel with the threshold value channel
    stage3_input = stage1_out.h5ad.combine(threshold_value_ch)
    
    // Unpack fields for process: tuple(study, sample, raw_h5ad), global_thresh
    processed_samples = ATAC_STAGE3_PROCESS(
        stage3_input.map { study, sample, h5ad, thresh -> tuple(study, sample, h5ad) },
        stage3_input.map { study, sample, h5ad, thresh -> thresh }
    )

    // 4. Final step: Group clean samples back together by study and generate the integrated atlas
    study_grouped_clean = processed_samples
        .groupTuple() // groups by key (study) automatically

    MERGE_AND_HARMONY_ATAC(study_grouped_clean)
}