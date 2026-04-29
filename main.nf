#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC_RAW {
    tag "${meta.id}"
    publishDir "${params.outdir}/qc/raw/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    """
}

process FASTP_TRIM {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastq_trimmed", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_trimmed_R1.fastq.gz"), path("${meta.id}_trimmed_R2.fastq.gz"), emit: trimmed
    path("${meta.id}.fastp.html"), emit: html
    path("${meta.id}.fastp.json"), emit: json

    script:
    """
    fastp \
      --in1 ${r1} \
      --in2 ${r2} \
      --out1 ${meta.id}_trimmed_R1.fastq.gz \
      --out2 ${meta.id}_trimmed_R2.fastq.gz \
      --html ${meta.id}.fastp.html \
      --json ${meta.id}.fastp.json \
      --thread ${task.cpus}
    """
}

process FASTQC_TRIMMED {
    tag "${meta.id}"
    publishDir "${params.outdir}/qc/trimmed/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    """
}

process BWA_MEM_MAP {
    tag "${meta.id}"
    publishDir "${params.outdir}/mapping/bam", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.mapped.bam"), path("${meta.id}.mapped.bam.bai"), path("${meta.id}.flagstat.tsv"), emit: bam

    script:
    def ref = params.bwa_index_prefix ?: params.reference_fasta
    """
    set -euo pipefail

    bwa mem -t ${task.cpus} ${ref} ${r1} ${r2} \
      | samtools collate -@ ${task.cpus} -O -u - \
      | samtools fixmate -@ ${task.cpus} -m -u - - \
      | samtools sort -@ ${task.cpus} -u - \
      | samtools markdup -@ ${task.cpus} - - \
      | samtools view -@ ${task.cpus} -b - \
      | samtools sort -@ ${task.cpus} -o ${meta.id}.mapped.bam -

    samtools index -@ ${task.cpus} ${meta.id}.mapped.bam
    samtools flagstat -@ ${task.cpus} -O tsv ${meta.id}.mapped.bam > ${meta.id}.flagstat.tsv
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_inputs

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}

process CNVKIT_BATCH {
    publishDir "${params.outdir}/cnvkit", mode: 'copy'

    input:
    path(case_bams)
    path(normal_bams)

    output:
    path("*.cnr"), emit: cnr
    path("*.cns"), emit: cns
    path("*.called.cns"), emit: called

    script:
    def seqMethod = (params.cnvkit_seq_method ?: 'wgs').toLowerCase()
    def caseArgs = case_bams.collect { it.toString() }.join(' ')
    def normalArgs = normal_bams.collect { it.toString() }.join(' ')
    def annotateArg = params.cnvkit_annotate ? "--annotate ${params.cnvkit_annotate}" : ''
    """
    cnvkit.py batch ${caseArgs} \
      --normal ${normalArgs} \
      --method ${seqMethod} \
      --fasta ${params.reference_fasta} \
      ${annotateArg} \
      --output-dir . \
      --processes ${task.cpus} \
      ${params.cnvkit_extra_batch_args}

    for cns in *.cns; do
      base=\$(basename "\${cns}" .cns)
      cnvkit.py call "\${cns}" --method ${params.cnvkit_call_method} --output "\${base}.called.cns"
    done
    """
}

workflow {
    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [id: row.sample, role: (row.cnv_role ?: '').toLowerCase()]
            tuple(meta, file(row.fastq_r1), file(row.fastq_r2))
        }

    FASTQC_RAW(samples_ch)
    FASTP_TRIM(samples_ch)
    FASTQC_TRIMMED(FASTP_TRIM.out.trimmed)
    BWA_MEM_MAP(FASTP_TRIM.out.trimmed)

    multiqc_inputs = FASTQC_RAW.out.html
        .mix(FASTQC_RAW.out.zip)
        .mix(FASTQC_TRIMMED.out.html)
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(FASTP_TRIM.out.html)
        .mix(FASTP_TRIM.out.json)
        .mix(BWA_MEM_MAP.out.bam.map { meta, bam, bai, flagstat -> flagstat })
        .collect()
    MULTIQC(multiqc_inputs)

    case_bams_ch = BWA_MEM_MAP.out.bam
        .filter { meta, bam, bai, flagstat -> !(meta.role in ['normal', 'control', 'reference']) }
        .map { meta, bam, bai, flagstat -> bam }
        .collect()

    normal_bams_ch = BWA_MEM_MAP.out.bam
        .filter { meta, bam, bai, flagstat -> meta.role in ['normal', 'control', 'reference'] }
        .map { meta, bam, bai, flagstat -> bam }
        .collect()

    CNVKIT_BATCH(case_bams_ch, normal_bams_ch)
}
