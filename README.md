# multiome-wnn-nf

A Nextflow (DSL2) pipeline for paired **snRNA-seq + snATAC-seq** multiome data:
SnapATAC2 (ATAC) and scanpy (RNA) feed into muon **WNN**, which is used to call
**broad cell types** (epithelial, fibroblast, endothelial, immune subtypes, …)
in colorectal cancer and healthy gut samples.

It does automatic per-sample **TSSe thresholding** (finds the dip), automatic
**marker-based annotation** from a config you control, and produces
**publication-quality plots that adapt to the number of samples** so text never
overlaps. Built to run on the KU Leuven **VSC (wice)** SLURM cluster.

---

## Pipeline

```
samplesheet ─┬─ ATAC_IMPORT (per sample)  fragments → tiles + TSSe + QC
             │        └─ ATAC_TSSE_THRESHOLD  auto "dip" per sample + global, QC plots
             │              └─ ATAC_EMBED  filter → spectral → harmony (+ gene activity)
             │
             └─ RNA_IMPORT (per sample)   10x .h5 → standardised AnnData
                      └─ RNA_QC_EMBED  concat → QC → normalise → PCA → harmony
                                              │
                              WNN_INTEGRATE  match barcodes → MuData → WNN → UMAP
                                              └─ ANNOTATE  leiden + marker scoring → broad_celltype
                                                       └─ PLOTS  UMAPs, dotplot, composition
```

Each box is a Nextflow process (`modules/`) wrapping a self-contained Python
script (`bin/`). Shared logic (adaptive plotting, the TSSe dip finder, a
shape-safe Harmony wrapper) lives in `bin/mwnn_utils.py`.

---

## Input

A CSV samplesheet with one row per multiome sample:

```csv
study,sample,atac_fragments,rna_h5
efremova,CRC01,/…/efremova/atac/efremova_CRC01_fragments.tsv.gz,/…/efremova/rna/efremova_CRC01_rna.h5
```

- `study` + `sample` → the internal id `<study>_<sample>` (e.g. `efremova_CRC01`),
  matching your `<study>_<sample>_…` naming scheme.
- The tabix index `<atac_fragments>.tbi` must sit next to the fragments file
  (it is staged automatically).
- Cell ids are standardised to `<study>_<sample>:<barcode>` on **both**
  modalities so RNA and ATAC line up cell-for-cell. Trailing `-1`-style
  barcode suffixes are stripped on both sides (`--strip_barcode_suffix`, on by
  default) to avoid the classic multiome matching mismatch.

---

## Running

```bash
# on VSC (wice), from a login node
module load Nextflow                 # or use your own nextflow binary
export NXF_WORK=$VSC_SCRATCH/nf-work # keep work dir on scratch

nextflow run main.nf \
    -profile vsc \
    --samplesheet assets/samplesheet.csv \
    --slurm_account <your_credit_account> \
    --outdir $VSC_SCRATCH/crc_multiome_results \
    --genome hg38 --species human
```

Quick local smoke test of wiring (no cluster):

```bash
nextflow run main.nf -profile test,conda --outdir results_test
```

Conda envs are activated per process. Defaults point at your existing kernels
(`snapatac2_jan26`, `scrna_nov25`, `muon_apr26`); override with
`--atac_env / --rna_env / --muon_env` (names or absolute env prefixes).

---

## The automatic TSSe threshold ("the dip")

`ATAC_TSSE_THRESHOLD` pools per-cell TSSe from every sample and finds the
**antimode** — the valley between the low-quality mode (dead/empty nuclei near
TSSe ≈ 1–2) and the high-quality mode — by KDE-smoothing each distribution,
locating its peaks and valleys, and taking the deepest valley between the
dominant low and high peaks. This is done **per sample** by default (quality
varies between samples) with a **global** dip computed too; set
`--tsse_mode global` to apply one threshold everywhere, or
`--tsse_mode manual --tsse_manual_threshold <value>` to skip auto-detection
and apply a fixed threshold everywhere (per-sample/global dips are still
computed and shown in `tsse_thresholds.csv`/the diagnostic plots for
reference, they just aren't the ones applied).

Guard rails, all configurable:
- `--tsse_min_high_peak` (4.0): the high-quality mode must sit above this.
- `--tsse_clamp_lo/hi` (2.0–10.0): thresholds are clamped into a safe band.
- `--default_min_tsse` (5.0): fallback when a distribution isn't clearly
  bimodal — and the `tsse_thresholds.csv` marks whether each sample was
  auto-detected or fell back, so nothing is silent.

Diagnostics written to `results/atac/qc/`: pooled histogram, per-sample KDE
ridge (▼ = detected dip), per-sample boxplot with the applied threshold, and a
TSSe-vs-fragments knee grid.

---

## RNA QC modes

`RNA_QC_EMBED` supports two filtering modes via `--rna_qc_mode`:

- **`basic`** (default): `--rna_min_genes`/`--rna_min_cells` filtering, plus
  an optional hard `--rna_max_pct_mt` cap (off by default — keeps all cells
  regardless of mito %).
- **`strict`**: everything in `basic`, plus the
  [sc-best-practices](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)
  outlier recipe — cells more than `--rna_strict_n_mads` (5) MADs from the
  median on log-total-counts, log-n-genes, or %-counts-in-top-20-genes are
  dropped, and mito is filtered by both a MAD threshold
  (`--rna_strict_mt_n_mads`, 3) and a hard cap (`--rna_strict_max_pct_mt`,
  8%). Both filters are logged with how many cells each removed.

---

## Automatic annotation

`ANNOTATE` clusters the WNN graph (Leiden, multiple resolutions), scores each
marker set from `assets/markers.yaml` on the RNA modality, and labels each
**cluster** by the highest mean z-scored marker score. A **margin rule**
(`params.margin`) makes ambiguous clusters `unknown` rather than a coin-flip,
and the ATAC gene-activity matrix gives an independent cross-check
(`atac_cluster_scores.csv`).

Everything is driven by `markers.yaml` — add a lineage or region set there, no
code changes:

```yaml
celltypes:
  fibroblast:
    rna: [DCN, COL1A1, COL1A2, COL3A1, PDGFRA, PDGFRB, LUM]
    atac_genes: [COL1A1, PDGFRA, DCN]
```

Outputs: `mdata_annotated.h5mu`, `cluster_assignments.csv`,
`cell_annotations.tsv` (per cell), and the score matrices.

---

## Adaptive plots

`bin/mwnn_utils.py` sizes figures, grids, palettes, point sizes, tick rotation
and legends by the number of samples / cells / cell types, so a 3-sample run
and a 40-sample run both come out readable. UMAP point size shrinks with cell
count; per-sample axes widen and rotate labels; categorical palettes extend
(glasbey when `colorcet` is installed). Final figures land in `results/plots/`.

---

## Key design decisions & assumptions

- **ATAC embedding = spectral-on-tiles** (SnapATAC2 native), not LSI-on-peaks.
  It needs no prior peak calling and is sufficient for broad WNN typing. Add a
  MACS3 peak-calling module downstream if you want peak-level analysis.
- **Harmony** is applied to both modalities over `sample`. The wrapper always
  returns a `(cells × dims)` array, fixing the `Z_corr` shape gotcha where the
  raw harmonypy output is `(dims × cells)`.
- **Broad typing is CPU-bound** (spectral, Harmony, scrublet, Leiden), so the
  `vsc` profile requests CPU nodes only. A commented `process_gpu` stub is in
  `conf/vsc.config` for when you add a GPU step (e.g. scVI) with its own H100
  credit account.
- Memory scales on retry (`conf/base.config`) so one oversized sample doesn't
  kill the run.

## What you'll likely tune first

- `assets/markers.yaml` — lineages/markers for your tissue.
- `--primary_resolution` and `--leiden_resolutions` — annotation granularity.
- `--tsse_clamp_lo/hi`, `--tsse_min_high_peak` — if a cohort's dip sits oddly.
- `conf/base.config` ceilings — to match wice node sizes and your credits.
