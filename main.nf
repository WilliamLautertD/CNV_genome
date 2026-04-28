#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def normalRoles = ['normal', 'control', 'reference'] as Set
def cnvMethods = params.cnv_method.tokenize(',').collect { it.trim().toLowerCase() } as Set
if (cnvMethods.contains('both')) {
    cnvMethods = ['cnvkit', 'gatk'] as Set
}
if (!cnvMethods.every { it in ['cnvkit', 'gatk'] }) {
    error "params.cnv_method must be 'cnvkit', 'gatk', or 'both'"
}

process FASTQC_RAW {
    tag "${meta.id}"
    publishDir "${params.outdir}/qc/raw/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.raw_fastqc.done")
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    touch ${meta.id}.raw_fastqc.done
    """

    stub:
    """
    touch ${meta.id}.raw_fastqc.done
    touch ${meta.id}_R1_fastqc.html ${meta.id}_R1_fastqc.zip
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

    stub:
    """
    touch ${meta.id}_trimmed_R1.fastq.gz
    touch ${meta.id}_trimmed_R2.fastq.gz
    touch ${meta.id}.fastp.html
    touch ${meta.id}.fastp.json
    """
}

process FASTQC_TRIMMED {
    tag "${meta.id}"
    publishDir "${params.outdir}/qc/trimmed/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.trimmed_fastqc.done")
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    touch ${meta.id}.trimmed_fastqc.done
    """

    stub:
    """
    touch ${meta.id}.trimmed_fastqc.done
    touch ${meta.id}_trimmed_R1_fastqc.html ${meta.id}_trimmed_R1_fastqc.zip
    """
}

process BWA_MEM_FILTER {
    tag "${meta.id}"
    publishDir "${params.outdir}/mapping/bam", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.marked.filtered.bam"), path("${meta.id}.marked.filtered.bam.bai"), path("${meta.id}.marked.filtered.flagstat.tsv"), emit: bam

    script:
    def ref = params.bwa_index_prefix ?: params.reference_fasta
    def readGroup = "@RG\\tID:${meta.id}\\tLB:${meta.library}\\tPL:${params.read_group_platform}\\tPU:${meta.unit}\\tSM:${meta.id}"
    """
    bwa mem -t ${task.cpus} -R '${readGroup}' ${ref} ${r1} ${r2} \
      | samtools collate -@ ${task.cpus} -O -u - \
      | samtools fixmate -@ ${task.cpus} -m -u - - \
      | samtools sort -@ ${task.cpus} -u - \
      | samtools markdup -@ ${task.cpus} - - \
      | samtools view -@ ${task.cpus} -b -q ${params.min_mapq} -F ${params.exclude_flags} - \
      | samtools sort -@ ${task.cpus} -o ${meta.id}.marked.filtered.bam -

    samtools index -@ ${task.cpus} ${meta.id}.marked.filtered.bam
    samtools flagstat -@ ${task.cpus} -O tsv ${meta.id}.marked.filtered.bam > ${meta.id}.marked.filtered.flagstat.tsv
    """

    stub:
    """
    touch ${meta.id}.marked.filtered.bam
    touch ${meta.id}.marked.filtered.bam.bai
    touch ${meta.id}.marked.filtered.flagstat.tsv
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

    stub:
    """
    touch multiqc_report.html
    """
}

process CNVKIT_REFERENCE {
    tag "${group}"
    publishDir "${params.outdir}/cnvkit/references", mode: 'copy'

    input:
    tuple val(group), path(normal_bams)

    output:
    tuple val(group), path("${group}.cnn"), emit: reference

    script:
    """
    cnvkit.py batch \
      --normal ${normal_bams} \
      --targets ${params.targets_bed} \
      --antitargets ${params.antitargets_bed} \
      --fasta ${params.reference_fasta} \
      --output-reference ${group}.cnn \
      --output-dir reference_work_${group} \
      --processes ${task.cpus}
    """

    stub:
    """
    touch ${group}.cnn
    """
}

process CNVKIT_BATCH {
    tag "${meta.id}"
    publishDir "${params.outdir}/cnvkit/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(flagstat)
    tuple val(group), path(reference)

    output:
    tuple val(meta), path("${meta.id}.marked.filtered.cns"), emit: cns

    script:
    def drop = params.cnvkit_drop_low_coverage ? '--drop-low-coverage' : ''
    """
    cnvkit.py batch ${bam} \
      --reference ${reference} \
      --output-dir . \
      --diagram \
      --scatter \
      --processes ${task.cpus} \
      ${drop} \
      ${params.cnvkit_extra_batch_args}
    """

    stub:
    """
    touch ${meta.id}.marked.filtered.cns
    """
}

process CNVKIT_CALL {
    tag "${meta.id}"
    publishDir "${params.outdir}/cnvkit/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(cns)

    output:
    tuple val(meta), path("${meta.id}.called.cns"), emit: called

    script:
    """
    cnvkit.py call ${cns} --method ${params.cnvkit_call_method} --output ${meta.id}.called.cns
    """

    stub:
    """
    touch ${meta.id}.called.cns
    """
}

process GATK_PREPROCESS_INTERVALS {
    publishDir "${params.outdir}/gatk_cnv/targets", mode: 'copy'

    output:
    tuple path("preprocessed_intervals.interval_list"), path("annotated_intervals.tsv"), emit: intervals

    script:
    """
    gatk PreprocessIntervals \
      -R ${params.reference_fasta} \
      -L ${params.targets_bed} \
      --interval-merging-rule ${params.gatk_interval_merging_rule} \
      -O preprocessed_intervals.interval_list

    gatk AnnotateIntervals \
      -R ${params.reference_fasta} \
      -L preprocessed_intervals.interval_list \
      --interval-merging-rule ${params.gatk_interval_merging_rule} \
      -O annotated_intervals.tsv
    """

    stub:
    """
    touch preprocessed_intervals.interval_list
    touch annotated_intervals.tsv
    """
}

process GATK_COLLECT_COUNTS {
    tag "${meta.id}"
    publishDir "${params.outdir}/gatk_cnv/counts", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(flagstat), path(intervals), path(annotated)

    output:
    tuple val(meta), path("${meta.id}.counts.hdf5"), emit: counts

    script:
    """
    gatk CollectReadCounts \
      -I ${bam} \
      -L ${intervals} \
      -R ${params.reference_fasta} \
      --interval-merging-rule ${params.gatk_interval_merging_rule} \
      --format HDF5 \
      -O ${meta.id}.counts.hdf5
    """

    stub:
    """
    touch ${meta.id}.counts.hdf5
    """
}

process GATK_PANEL_OF_NORMALS {
    tag "${group}"
    publishDir "${params.outdir}/gatk_cnv/pon", mode: 'copy'

    input:
    tuple val(group), path(normal_counts), path(annotated)

    output:
    tuple val(group), path("${group}.pon.hdf5"), emit: pon

    script:
    def countArgs = normal_counts.collect { "-I ${it}" }.join(' ')
    """
    gatk CreateReadCountPanelOfNormals \
      ${countArgs} \
      --annotated-intervals ${annotated} \
      -O ${group}.pon.hdf5
    """

    stub:
    """
    touch ${group}.pon.hdf5
    """
}

process GATK_DENOISE {
    tag "${meta.id}"
    publishDir "${params.outdir}/gatk_cnv/denoised", mode: 'copy'

    input:
    tuple val(meta), path(counts)
    tuple val(group), path(pon)

    output:
    tuple val(meta), path("${meta.id}.standardizedCR.tsv"), path("${meta.id}.denoisedCR.tsv"), emit: denoised

    script:
    """
    gatk DenoiseReadCounts \
      -I ${counts} \
      --count-panel-of-normals ${pon} \
      --standardized-copy-ratios ${meta.id}.standardizedCR.tsv \
      --denoised-copy-ratios ${meta.id}.denoisedCR.tsv
    """

    stub:
    """
    touch ${meta.id}.standardizedCR.tsv
    touch ${meta.id}.denoisedCR.tsv
    """
}

process GATK_MODEL_SEGMENTS {
    tag "${meta.id}"
    publishDir "${params.outdir}/gatk_cnv/segments", mode: 'copy'

    input:
    tuple val(meta), path(standardized), path(denoised)

    output:
    tuple val(meta), path("${meta.id}.cr.seg"), emit: segments

    script:
    """
    gatk ModelSegments \
      --denoised-copy-ratios ${denoised} \
      --output . \
      --output-prefix ${meta.id} \
      --number-of-changepoints-penalty-factor ${params.gatk_changepoint_penalty}
    """

    stub:
    """
    touch ${meta.id}.cr.seg
    """
}

process GATK_CALL_SEGMENTS {
    tag "${meta.id}"
    publishDir "${params.outdir}/gatk_cnv/segments", mode: 'copy'

    input:
    tuple val(meta), path(segments)

    output:
    tuple val(meta), path("${meta.id}.called.seg"), emit: called

    script:
    """
    gatk CallCopyRatioSegments --input ${segments} --output ${meta.id}.called.seg
    """

    stub:
    """
    touch ${meta.id}.called.seg
    """
}

workflow {
    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id: row.sample,
                role: row.cnv_role,
                group: row.cnv_reference_group,
                library: row.read_group_library ?: 'lib1',
                unit: row.read_group_unit ?: 'unit1'
            ]
            tuple(meta, file(row.fastq_r1), file(row.fastq_r2))
        }

    raw_reads_ch = samples_ch.map { meta, r1, r2 -> tuple(meta, [r1, r2]) }

    FASTQC_RAW(raw_reads_ch)
    FASTP_TRIM(samples_ch)
    FASTQC_TRIMMED(FASTP_TRIM.out.trimmed)
    BWA_MEM_FILTER(FASTP_TRIM.out.trimmed)

    multiqc_inputs = FASTQC_RAW.out.html
        .mix(FASTQC_RAW.out.zip)
        .mix(FASTQC_TRIMMED.out.html)
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(FASTP_TRIM.out.html)
        .mix(FASTP_TRIM.out.json)
        .mix(BWA_MEM_FILTER.out.bam.map { meta, bam, bai, flagstat -> flagstat })
        .collect()
    MULTIQC(multiqc_inputs)

    if (cnvMethods.contains('cnvkit')) {
        normal_bams_by_group = BWA_MEM_FILTER.out.bam
            .filter { meta, bam, bai, flagstat -> normalRoles.contains(meta.role.toLowerCase()) }
            .map { meta, bam, bai, flagstat -> tuple(meta.group, bam) }
            .groupTuple()

        CNVKIT_REFERENCE(normal_bams_by_group)

        case_bams_by_group = BWA_MEM_FILTER.out.bam
            .filter { meta, bam, bai, flagstat -> !normalRoles.contains(meta.role.toLowerCase()) }
            .map { meta, bam, bai, flagstat -> tuple(meta.group, meta, bam, bai, flagstat) }

        cnvkit_refs_by_group = CNVKIT_REFERENCE.out.reference
            .map { group, reference -> tuple(group, reference) }

        case_bams_by_group
            .join(cnvkit_refs_by_group)
            .map { group, meta, bam, bai, flagstat, reference ->
                tuple(tuple(meta, bam, bai, flagstat), tuple(group, reference))
            }
            .multiMap { sample_tuple, ref_tuple ->
                sample: sample_tuple
                ref: ref_tuple
            }
            .set { cnvkit_inputs }

        CNVKIT_BATCH(cnvkit_inputs.sample, cnvkit_inputs.ref)
        CNVKIT_CALL(CNVKIT_BATCH.out.cns)
    }

    if (cnvMethods.contains('gatk')) {
        GATK_PREPROCESS_INTERVALS()
        gatk_count_inputs = BWA_MEM_FILTER.out.bam.combine(GATK_PREPROCESS_INTERVALS.out.intervals)
        GATK_COLLECT_COUNTS(gatk_count_inputs)

        normal_counts_by_group = GATK_COLLECT_COUNTS.out.counts
            .filter { meta, counts -> normalRoles.contains(meta.role.toLowerCase()) }
            .map { meta, counts -> tuple(meta.group, counts) }
            .groupTuple()

        annotated_intervals = GATK_PREPROCESS_INTERVALS.out.intervals.map { intervals, annotated -> annotated }
        normal_counts_with_intervals = normal_counts_by_group.combine(annotated_intervals)
        GATK_PANEL_OF_NORMALS(normal_counts_with_intervals)

        case_counts_by_group = GATK_COLLECT_COUNTS.out.counts
            .filter { meta, counts -> !normalRoles.contains(meta.role.toLowerCase()) }
            .map { meta, counts -> tuple(meta.group, meta, counts) }

        gatk_pons_by_group = GATK_PANEL_OF_NORMALS.out.pon
            .map { group, pon -> tuple(group, pon) }

        case_counts_by_group
            .join(gatk_pons_by_group)
            .map { group, meta, counts, pon ->
                tuple(tuple(meta, counts), tuple(group, pon))
            }
            .multiMap { counts_tuple, pon_tuple ->
                counts: counts_tuple
                pon: pon_tuple
            }
            .set { gatk_inputs }

        GATK_DENOISE(gatk_inputs.counts, gatk_inputs.pon)
        GATK_MODEL_SEGMENTS(GATK_DENOISE.out.denoised)
        GATK_CALL_SEGMENTS(GATK_MODEL_SEGMENTS.out.segments)
    }
}
