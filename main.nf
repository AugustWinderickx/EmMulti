#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ATAC_IMPORT; ATAC_TSSE_THRESHOLD; ATAC_EMBED } from './modules/atac.nf'
include { RNA_IMPORT; RNA_QC_EMBED }                     from './modules/rna.nf'
include { WNN_INTEGRATE; ANNOTATE; PLOTS }               from './modules/wnn.nf'

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
      --tsse_mode            per_sample | global   (auto TSSe dip)
      --markers              marker YAML (default assets/markers.yaml)
      --primary_resolution   leiden resolution used for annotation
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
    RNA_IMPORT(rna_in)
    RNA_QC_EMBED(RNA_IMPORT.out.h5ad.collect())

    // ---- WNN + annotation + plots ----------------------------------------
    WNN_INTEGRATE(RNA_QC_EMBED.out.rna, ATAC_EMBED.out.atac)
    

    markers_ch = Channel.value(file(params.markers, checkIfExists: true))
    ANNOTATE(WNN_INTEGRATE.out.mdata, markers_ch, ATAC_EMBED.out.gene_activity)
    PLOTS(ANNOTATE.out.mdata, markers_ch, ATAC_EMBED.out.gene_activity)
}

workflow.onComplete {
    log.info (workflow.success
        ? "\n✅  Done. Results in ${params.outdir}\n"
        : "\n❌  Failed after ${workflow.duration}\n")
}
