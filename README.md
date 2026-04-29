# CNV_genome

Simple Nextflow pipeline for:
1. FastQC on raw reads
2. fastp trimming
3. FastQC on trimmed reads
4. BWA mapping + BAM sorting/indexing
5. MultiQC summary
6. CNVkit batch run

## Required input

Edit `config/samples.tsv` with one sample per line:
- `sample`
- `fastq_r1`
- `fastq_r2`
- `cnv_role` (`normal`/`control`/`reference` for normals; anything else is treated as case)

## Configure

Edit `nextflow.config`:
- `params.reference_fasta`
- optional `params.bwa_index_prefix`
- `params.cnvkit_seq_method` (default `wgs`)
- optional `params.cnvkit_annotate` for gene annotation (refFlat)

## Run

```bash
conda env create -f env/qc_mapping_cnv.yaml
conda activate qc_mapping_cnv
nextflow run main.nf -profile conda -resume
```

## Key outputs

- Mapped BAMs: `results/automated_pipeline/mapping/bam/*.mapped.bam`
- BAM indexes: `results/automated_pipeline/mapping/bam/*.mapped.bam.bai`
- Mapping stats: `results/automated_pipeline/mapping/bam/*.flagstat.tsv`
- CNVkit outputs: `results/automated_pipeline/cnvkit/*.cnr`, `*.cns`, and `*.called.cns`
- MultiQC report: `results/automated_pipeline/qc/multiqc_report.html`
