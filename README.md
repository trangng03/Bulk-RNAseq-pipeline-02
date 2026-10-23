# Bulk RNA-seq: STAR + kallisto + RSeQC + DESeq2 + GSEA + rMATS

For paired-end bulk RNA-seq, two conditions and at least two independent biological replicates per condition. Example sheet: control1–3 and NORADkd1–3 (Condition labels `control` and `treated`). Sample IDs and condition labels are different: `case` must match the Condition column, not an individual sample ID.

Nothing installs packages, creates environments, downloads data, or activates environments automatically. Use the existing `rnaseq` environment for QC, kallisto and R analysis; use the existing `rmats` environment for rMATS if needed. This is a new package; it does not modify your earlier pipeline.

## Workflow

FASTQs → FastQC → STAR → RSeQC QC and strandedness assessment.
The same FASTQs → kallisto (using assessed strandedness) → tximport → DESeq2 → GSEA.
STAR genome BAMs → rMATS, independently of the expression analysis.
MultiQC summarizes supported QC outputs. rMATS results are summarized separately.

STAR GeneCounts are saved for reference but DESeq2 in this package uses kallisto estimates imported through tximport. TPM is not supplied directly to DESeq2. `DESeqDataSetFromTximport` handles estimated counts and transcript-length offsets.

## Files in this package

- `run.sh`: Bash entry point for each stage.
- `config.json`: paths, comparison and settings to edit.
- `samples.tsv`: six-sample template to edit; actual tab-separated fields.
- `check_dependencies.sh`: inspect commands and R packages without installing.
- `build_kallisto_index.sh`: optional index builder from your existing transcript FASTA.
- `scripts/pipeline.py`: standard-library Python helper for annotation preparation, commands, QC, checkpoints and rMATS summaries.
- `scripts/analyze.R`: tximport, DESeq2, gene-name mapping, expression plots and fgsea.
- `scripts/test_helpers.py`: small fixture tests for annotation conversion, strand interpretation, commands and event filtering.

## 1. Set your paths

Extract the ZIP and open a terminal in `rnaseq_pipeline`.
Edit `config.json`, replacing all `/EDIT/...` entries. Relative paths are resolved relative to the config file, including FASTQ paths from the sample sheet. Use a NEW results directory for this pipeline, such as `/home/tnguyen25/Documents/Projects/RNAseq/results_kallisto_rmats`, so it does not overwrite your earlier results.

Required references:

1. STAR genome index directory containing `genomeParameters.txt` (only needed when aligning anew).
2. kallisto transcript index FILE, not a directory.
3. Annotation GTF, plain or gzipped. It must include exon records with transcript_id and gene_id.
4. Human gene-symbol GMT file, for example your Hallmark GMT.

Use genome and transcript sequences, GTF and indexes from the same reference release/build, with matching chromosome names. The package derives BED12 transcript models (including coding boundaries when CDS records exist), transcript-to-gene mapping, and gene-symbol lookup from that GTF. It creates a plain GTF copy for rMATS. It does not verify the biological identity of prebuilt indexes; verify their reference provenance yourself.

The kallisto transcript FASTA must use transcript IDs compatible with the GTF. Ensembl transcript IDs with optional terminal numeric versions are supported. GENCODE headers containing pipes need preprocessing or corresponding changes to tximport; do not mix incompatible headers and mappings. Include the transcript types relevant to your experiment, especially noncoding transcripts for a NORAD study.

Optional index creation, if you do not already have a compatible index:

```bash
conda activate rnaseq
bash build_kallisto_index.sh /absolute/path/transcripts.fa.gz /absolute/path/new_transcripts.idx
```

The script refuses to overwrite an existing index. Build your STAR genome index separately if you do not have one.

Edit `samples.tsv` to contain the actual R1/R2 paths, one row per biological replicate:

```text
Sample  Condition  Read1  Read2
control1  control  /path/control1_1.fastq.gz  /path/control1_2.fastq.gz
NORADkd1  treated  /path/NORADkd1_1.fastq.gz  /path/NORADkd1_2.fastq.gz
```

Use TABS, not spaces, in the actual sheet. Do not use the illustrative two-row example as the complete dataset. FASTQs may be plain or gzip-compressed; mates must have the same compression type. R1 and R2 must not be interchanged. Keep all lanes for a biological replicate together before using this one-pair-per-sample template.

Configuration settings:

- `threads`: default 8, adjust to your CPU/RAM capacity. Samples run sequentially.
- `case`: `treated`; `control`: `control` in the provided sheet.
- `strand`: `auto` by default; alternatives `reverse`, `forward`, `unstranded`.
- `read_length`: 99 in the template based on your earlier dataset; verify it for your actual reads. QC records lengths in the first 1,000 reads of each mate in `03_rseqc/read_lengths.tsv`. Unless variable_read_length is enabled, rMATS preparation stops when these sampled lengths differ from the configured length. This sample check cannot prove that every read has that length.
- `variable_read_length`: false; use true for variable-length reads, keeping read_length representative of your sequencing reads. See rMATS documentation.
- `allow_clipping`: false; rMATS excludes clipped alignments by default. Change deliberately after considering the rMATS documentation.
- `rmats_fdr`: 0.05; `rmats_delta_psi`: 0.1 (10 percentage points).
- `rmats_command`: a JSON array, e.g. `["rmats.py"]`, `["python", "/absolute/path/rmats.py"]`, or `["/absolute/path/run_rmats"]`. Do not put the entire command in one string.

## 2. Check your current environment

```bash
conda activate rnaseq
bash check_dependencies.sh rnaseq
```

Required commands: Python3, Rscript, STAR, samtools, FastQC, MultiQC, kallisto, infer_experiment.py, read_distribution.py, geneBody_coverage.py and junction_annotation.py. STAR is unnecessary if reusing existing BAMs; the broad dependency check still reports it. RSeQC gene-body/junction plotting also requires working R in the active environment.

Required R packages: tximport, DESeq2, fgsea, ggplot2, ggrepel, SummarizedExperiment.

If something is missing, the scripts report it and stop. Nothing is installed. A package available in RStudio might not be available to the terminal's Rscript.

## 3. Run QC and inspect results

```bash
bash run.sh qc config.json
```

This runs FastQC, STAR (unless reusing BAMs), samtools BAM checks/indexing/flagstat and four RSeQC checks per sample. It also generates a MultiQC report.

Inspect:

- `results/08_multiqc/multiqc_report.html`
- `results/03_rseqc/strandedness.tsv`
- Per-sample `infer_experiment.txt`, read distribution, gene-body coverage and junction annotation outputs.

The strand heuristic normalizes the two orientation fractions to their sum. It requires at least 50% of sampled reads to be explainable. A forward proportion >=0.8 suggests forward; <=0.2 suggests reverse; 0.4–0.6 suggests unstranded; other cases are unclear. These are conservative workflow rules, not official universal cutoffs. `auto` continues only if ALL samples have the same clear suggestion. It stops before kallisto/rMATS when suggestions are unclear or inconsistent. Review the library protocol and individual results before setting an explicit strand. An explicit strand applies to every sample; mixed library protocols need separate runs.

For your previously measured samples, reverse is likely, but the actual RSeQC report is what this package checks.

Strand translation:

| Setting | kallisto | rMATS |
|---|---|---|
| unstranded | no strand flag | fr-unstranded |
| forward | --fr-stranded | fr-secondstrand |
| reverse | --rf-stranded | fr-firststrand |

These names differ from the numeric STAR GeneCounts convention.

### Reuse your existing STAR BAMs

Set this in config.json BEFORE running QC:

```json
"existing_star_dir": "/home/tnguyen25/Documents/Projects/RNAseq/results/02_star"
```

Each sample must be under `<existing_star_dir>/<Sample>/Aligned.sortedByCoord.out.bam`. The package checks and indexes those BAMs in place, but does not realign them. They must match the sample sheet and reference, and contain splice-aware genome alignments. Reuse requires the BAM directory to be writable for indexing.

## 4. Run kallisto and expression analysis

```bash
bash run.sh quant config.json
bash run.sh analysis config.json
```

The R script reads kallisto abundance.tsv files and the GTF-derived transcript-to-gene table. It retains genes with >=10 counts in at least two samples. DESeq2 compares case/control. Significant-gene CSV uses padj<0.05 and |log2FoldChange|>=1. The top up/down volcano labels select the largest positive and most negative fold changes among qualifying genes, not the smallest p-values. Fewer than five qualifying genes in a direction yields fewer labels.

CSV results use gene symbols, falling back to Ensembl IDs if symbols are absent. Original IDs remain available in `gene_name_mapping.csv` and the RDS object. Different IDs can share a symbol; expression CSVs do not collapse them. GSEA requires one statistic per symbol, so it keeps the largest absolute Wald statistic when symbols repeat. GSEA uses all finite mapped statistics from retained genes, not only significant genes.

Expression outputs: all/significant genes, normalized counts, DESeq2 object, MA, PCA, sample correlation, top-30 significant-gene heatmap, volcano PNG/PDF and labeled-gene CSV. The heatmap is omitted if fewer than two qualifying variable genes exist. Sample PCA/correlations use variance-stabilized expression, not TPM.

GSEA outputs: ranked genes, all/significant enrichment results, a top-pathway NES bar plot and up to six running-enrichment plots. Top-pathway plots rank by padj and can include nonsignificant pathways; inspect the table. Positive NES corresponds to the positive case/control ranking; negative NES to its negative end. Neither directly proves protein activity changes.

## 5. Run rMATS in your existing rMATS environment

```bash
conda activate rmats
bash check_dependencies.sh rmats
bash run.sh rmats config.json
```

Check your configured rMATS command with `--help` before the full run. If the command is `python /path/rmats.py`, use that form for the help check as well. Activation is manual; no environment is created.

Group 1 (`--b1`) is case; group 2 (`--b2`) is control. Positive IncLevelDifference means higher event inclusion in case; negative means lower inclusion. This is event-specific and is NOT a gene expression fold change. Paired-end sequencing is enabled, but this package uses independent biological groups, not rMATS' paired-subject statistical model.

rMATS produces SE, A5SS, A3SS, MXE and RI event results. Both JC (junction counts) and JCEC (junction+exon counts) remain available; do not add their event totals together as independent discoveries. Filtered tables require FDR<0.05 and |IncLevelDifference|>=0.1 by default. The package also writes a per-table event summary and records group direction.

For chosen events, sashimi plots/IGV inspection are useful follow-ups, but are not automated here. Trimming, batch adjustment, isoform-switch analysis and additional enrichment of splicing genes are also outside this package. Do not trim blindly based only on low-quality tails; if trimming is needed, update the FASTQ paths and read-length settings and use a fresh results directory.

## 6. Refresh the report

```bash
conda activate rnaseq
bash run.sh report config.json
```

You can run all stages with `bash run.sh all config.json` ONLY when all tools (including rMATS) are already available in that one active environment. Separate stage commands are recommended for your rnaseq/rmats setup. No GSEA output is used as rMATS input.

## Output folders

| Folder | Contents |
|---|---|
| 00_reference | Plain GTF copy, BED12, tx2gene and symbol lookup |
| 01_fastqc | Per-sample raw FASTQ reports |
| 02_star | New STAR alignments/GeneCounts; empty if external BAMs reused |
| 03_rseqc | Alignment QC and strand summary |
| 04_kallisto | Per-sample abundance, estimated gene counts, strand used |
| 05_deseq2 | Expression results and plots |
| 06_gsea | Enrichment results and plots |
| 07_rmats | Raw event files, filtered events, summary and group BAM lists |
| 08_multiqc | Combined QC HTML report |
| logs | Command logs, executable paths, R session info, stage completion records |

The root also contains input_manifest.tsv and resolved_config.json.

## Reruns and logs

Completed stages have completion markers. A repeated stage with the same fingerprint is skipped. Changed config, scripts or input file metadata causes an error rather than silent reuse. To deliberately rerun use `--force`, e.g. `bash run.sh analysis config.json --force`. Existing expression CSVs/plots are overwritten. This is stage-level resumption: an interrupted QC stage may realign earlier samples on rerun. Never delete result files and expect a completed marker to regenerate them without --force.

rMATS refuses to reuse a nonempty raw result directory. Move `07_rmats/raw` and `07_rmats/tmp` out of the results tree or use a new results_dir before deliberately rerunning it. `--force` does not bypass that check. No folders are deleted automatically.

On a tool error, inspect the log path printed in the terminal. Full biological/software compatibility depends on your installed versions and references.

## Validation supplied with this package

```bash
bash -n run.sh
bash -n check_dependencies.sh
python3 scripts/test_helpers.py
Rscript -e 'parse(file="scripts/analyze.R"); cat("R syntax OK\\n")'
```

During package creation: shell syntax checks and Python compilation passed. Three fixture tests passed, covering BED/transcript mapping, strandedness classification, stage argument wiring and rMATS event filtering. The bioinformatics commands were simulated for wiring tests, not run on real sequencing data. No R executable or real reference/FASTQ data were available in the creation environment, so the R script and end-to-end biological analysis remain to be executed in your environment. Do not treat the fixture tests as validation of biological conclusions.

## Official documentation

- kallisto: https://pachterlab.github.io/kallisto/manual
- tximport/DESeq2: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
- RSeQC: https://rseqc.sourceforge.net/
- rMATS: https://github.com/Xinglab/rmats-turbo
- fgsea: https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
