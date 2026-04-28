# CNV_genome
Repository for automated FASTQ to CNV analysis
## QC pipelines
deeptools tools for exploring deep sequencing data
- Analysis of correlation of bamfiles

## Mapping pipelines

## Automated FASTQ to CNV workflow

An end-to-end Nextflow workflow now lives in `main.nf`. It automates:

1. Raw-read FastQC
2. fastp trimming
3. Trimmed-read FastQC
4. BWA-MEM mapping with read groups
5. samtools duplicate marking and filtered BAM generation with `samtools flagstat`
6. MultiQC summary
7. CNV calling with CNVkit, GATK CNV, or both

### Files

- `config/samples.tsv` - one row per sample, with FASTQ paths and CNV reference group.
- `nextflow.config` - reference genome, targets, output directory, executors, resources, and CNV settings.
- `envs/qc_mapping_cnv.yaml` - Conda environment for Nextflow and the required tools.
- `main.nf` - automated pipeline definition.

### Quick start

Edit `config/samples.tsv` and `nextflow.config`, then run:

```bash
conda env create -f envs/qc_mapping_cnv.yaml
conda activate qc_mapping_cnv
nextflow run main.nf -profile conda -resume
```

For an HPC run using the built-in SLURM profile:

```bash
nextflow run main.nf -profile slurm -resume
```

### CNV modes

Set `params.cnv_method` in `nextflow.config`:

- `cnvkit` - build CNVkit references from samples marked `normal`, `control`, or `reference`, then call `.called.cns` files.
- `gatk` - preprocess target intervals, collect read counts, build a panel of normals, denoise, segment, and call `.called.seg` files.
- `both` - run both CNVkit and GATK CNV outputs from the same marked BAMs.

### CNVkit WGS support

The CNVkit path is now WGS-first (`params.cnvkit_seq_method = 'wgs'`) so it works for whole-genome data instead of only WES/hybrid capture assumptions.

- In WGS mode, the workflow runs `cnvkit.py batch --method wgs` and uses sequencing-accessible regions (`--access`) when provided; otherwise CNVkit computes access regions from the FASTA internally.
- In hybrid/amplicon mode, it uses `--targets` and `--antitargets` as before.
- `params.cnvkit_no_edge = true` is enabled by default, matching CNVkit WGS recommendations.
- To enable gene labels, either set `params.cnvkit_annotate` directly to the refFlat path, or set `params.cnvkit_annotate = 'true'` and provide `params.cnvkit_refflat = '/refs/refFlat.txt'`.
- Backward compatibility: `params.annotate` and `params.refflat` are accepted aliases.

Samples with `cnv_role` set to `normal`, `control`, or `reference` are used for CNV references. Other roles, such as `case` or `treated`, are CNV-called against their `cnv_reference_group`.

## Generic ChIP-seq Nextflow workflow

For ChIP-seq, a separate generic Nextflow workflow is available in `chipseq_main.nf` with configuration in `chipseq_nextflow.config`.

### Why this is generic

- Input FASTQ files are fully controlled by `config/chipseq_samples.tsv` (no fixed filename pattern is required).
- Sample names can be any value in the `sample` column.
- The `control_sample` column pairs each ChIP sample with the exact Input sample ID for downstream `bamCompare` output.
- Paths and run parameters are set in `chipseq_nextflow.config`, including reference genome, mapping filters, bigWig options, and `bamCompare` settings (ratio mode with duplicate ignoring, no scale-factor normalization).


## Copy Number Variation (CNV) Analysis Pipelines
GATK & CNVkit Workflows for Targeted and Whole-Exome Sequencing

### Overview
This repository provides reproducible, HPC-ready workflows for copy number variation (CNV) analysis using two independent pipelines:
1. GATK CNV Workflow — Best-practice CNV calling using the Broad Institute's Genome Analysis Toolkit (GATK).
2. CNVkit Workflow — Coverage-based CNV detection using CNVkit for targeted and hybrid capture sequencing.

Each workflow includes:
- Ready-to-run SLURM batch scripts for HPC clusters
- Step-by-step setup and execution guides
- Notes on parameters, expected outputs, and biological interpretation

### Reproducibility
All scripts are fully modular and can be customized per project.
Each step includes:
- Input and output definitions
- Environment setup instructions
- Optional parameters for advanced tuning
To rerun or adapt:
- Update paths in the scripts (BAM, REF, TARGETS, etc.)
- Submit each job to the HPC queue using sbatch
- Review logs and resulting CNV tables/plots

## Duplicate and fusion genes
- Manual inspection of fusioned genes 
- Using the "supplementary", "mates on different chromosomes", and mates on same chromosomes but in distant than expected" reads. 
- Compare it with Normal. 

## Repository Structure

### Folder Organization

- **CNVkit**: Contains scripts and tools for copy number variation analysis using CNVkit.
- **GATK_CNV**: Includes files related to the Genome Analysis Toolkit for copy number variations.
- **Mapping**: Houses the mapping files and scripts used for aligning sequencing data.
- **QC**: Contains quality control metrics and reports for the datasets.
- **ChIP-Seq_Chromatin_analysis**: Includes analysis scripts and data related to ChIP-Seq experiments.
- **Duplication_fusion_genes**: Contains files related to the analysis of gene duplications and fusions.
- **deeptools**: Houses scripts and tools used for deep data analysis.
- **bcftools**: Contains tools for variant calling and manipulating VCF files.
- **fastq**: Houses FASTQ files of raw sequencing data.


## Citation
relevant tools:
- GATK CNV – Benjamin et al., Nature Genetics (2013)
- CNVkit – Talevich et al., PLOS Computational Biology (2016)
