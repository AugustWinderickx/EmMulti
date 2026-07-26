#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ATAC_IMPORT; ATAC_TSSE_THRESHOLD; ATAC_EMBED } from './modules/atac.nf'
include { RNA_IMPORT; RNA_QC_EMBED }                     from './modules/rna.nf'
include { WNN_INTEGRATE as WNN_INTEGRATE_HARMONY
          WNN_INTEGRATE as WNN_INTEGRATE_RAW
          ANNOTATE      as ANNOTATE_HARMONY
          ANNOTATE      as ANNOTATE_RAW
          PLOTS         as PLOTS_HARMONY
          PLOTS         as PLOTS_RAW }                   from './modules/wnn.nf'

def helpMessage() {
    log.info """
    multiome-wnn-nf — snRNA + snATAC (SnapATAC2 -> scanpy -> muon/WNN) broad
    cell-type annotation for CRC / healthy gut.

    Usage:
      nextflow run main.nf -profile vsc --samplesheet assets/samplesheet.csv \\
          --outdir results --genome hg38

    Key params (see nextflow.config for all + defaults):
      --samplesheet          CSV: study,sample,atac_fragments,rna_h5
      --outdir               results directory
      --genome               hg38 | mm10
      --species              human | mouse
      --gtf                  GENCODE/Ensembl GTF used to fill in gene
                              symbols for RNA inputs that carry Ensembl IDs
                              as var_names instead of symbols
      --tsse_mode            per_sample | global | manual   (auto TSSe dip,
                              or manual with --tsse_manual_threshold)
      --rna_qc_mode          basic | strict   (strict adds MAD-based outlier
                              + mito filtering, see nextflow.config)
      --markers              marker YAML (default assets/markers.yaml)
      --primary_resolution   leiden resolution used for annotation
      --barcode_translation  auto | force | skip -- per-sample 10x Multiome
                              ATAC->GEX barcode translation when direct
                              matching is poor (auto also flags samples with
                              no true multiome pairing; see nextflow.config)
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); exit 0 }

    // ---- parse samplesheet -> per-sample channels -------------------------
    Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [study: row.study, sample: row.sample,
                        id: "${row.study}_${row.sample}"]
            def frag = file(row.atac_fragments, checkIfExists: true)
            def tbi  = file("${row.atac_fragments}.tbi")   // tabix index alongside
            def rna  = file(row.rna_h5, checkIfExists: true)
            [meta, frag, tbi, rna]
        }
        .set { samples }

    atac_in = samples.map { meta, frag, tbi, rna -> [meta, frag, tbi] }
    rna_in  = samples.map { meta, frag, tbi, rna -> [meta, rna] }

    // ---- ATAC branch ------------------------------------------------------
    ATAC_IMPORT(atac_in)
    ATAC_TSSE_THRESHOLD(ATAC_IMPORT.out.qc.collect())
    ATAC_EMBED(ATAC_IMPORT.out.h5ad.collect(), ATAC_TSSE_THRESHOLD.out.thresholds)

    // ---- RNA branch -------------------------------------------------------
    gtf_ch = Channel.value(file(params.gtf, checkIfExists: true))
    RNA_IMPORT(rna_in, gtf_ch)
    RNA_QC_EMBED(RNA_IMPORT.out.h5ad.collect())

    // ---- WNN + annotation + plots ------------------------------------------
    // Run the Harmony-corrected branch and, in parallel, an uncorrected
    // ("raw") branch on the same RNA/ATAC embeddings so batch-correction
    // effects can be compared side by side. Each process is aliased per
    // variant on import (DSL2 forbids invoking the same process twice in one
    // workflow), so the two branches run as independent process instances.
    markers_ch = Channel.value(file(params.markers, checkIfExists: true))
    atac_to_gex_translation_ch = Channel.value(
        file(params.atac_to_gex_translation, checkIfExists: true))

    WNN_INTEGRATE_HARMONY(RNA_QC_EMBED.out.rna, ATAC_EMBED.out.atac,
                           'harmony', params.rna_rep, params.atac_rep,
                           atac_to_gex_translation_ch)
    ANNOTATE_HARMONY(WNN_INTEGRATE_HARMONY.out.mdata, markers_ch,
                      ATAC_EMBED.out.gene_activity, 'harmony')
    PLOTS_HARMONY(ANNOTATE_HARMONY.out.mdata, markers_ch,
                   ATAC_EMBED.out.gene_activity, 'harmony')

    if (params.compare_batch_correction) {
        WNN_INTEGRATE_RAW(RNA_QC_EMBED.out.rna, ATAC_EMBED.out.atac,
                           'raw', params.rna_rep_raw, params.atac_rep_raw,
                           atac_to_gex_translation_ch)
        ANNOTATE_RAW(WNN_INTEGRATE_RAW.out.mdata, markers_ch,
                      ATAC_EMBED.out.gene_activity, 'raw')
        PLOTS_RAW(ANNOTATE_RAW.out.mdata, markers_ch,
                   ATAC_EMBED.out.gene_activity, 'raw')
    }
}

workflow.onComplete {
    log.info (workflow.success
        ? "\n✅  Done. Results in ${params.outdir}\n"
        : "\n❌  Failed after ${workflow.duration}\n")
}
