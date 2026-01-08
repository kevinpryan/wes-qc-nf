nextflow.enable.dsl = 2

/*
 * -------------------------------------------------
 *  WES QC PIPELINE (with IDT Exome v2 specific intervals)
 * -------------------------------------------------
 *  Metrics: Chimeras, Unmapped %, On-Target,
 *           Fold 80, FPR (via Ti/Tv)
 * -------------------------------------------------
 */

// Default parameters (can be overwritten via command line)
params.bams   = "data/*.bam"            // Path to input BAMs
params.ref    = "refs/genome.fasta"     // Path to Reference Fasta
params.probes_bed = "refs/probes.hg38.bed" // Path to IDT Probes BED file
params.targets_bed = "refs/targets.hg38.bed" // Path to IDT Targets BED file
params.outdir = "results"               // Output directory

log.info """
W E S   Q C   P I P E L I N E (IDT Exome v2)
============================================
BAMs        : ${params.bams}
Reference   : ${params.ref}
Probes BED  : ${params.probes_bed}
Targets BED : ${params.targets_bed}
Output      : ${params.outdir}
"""

process SAMTOOLS_SORT {
    label "SAMTOOLS_CONTAINER"
    tag "$bam_id"
    // We don't necessarily need to publish the sorted BAMs if you only want metrics, 
    // but it's good practice to keep them.
    publishDir "${params.outdir}/sorted_bams", mode: 'copy'

    input:
    tuple val(bam_id), path(bam)

    output:
    tuple val(bam_id), path("${bam_id}.sorted.bam"), emit: bam
    path "${bam_id}.sorted.bam.bai", emit: index

    script:
    """
    samtools sort -@ ${task.cpus} -o ${bam_id}.sorted.bam $bam
    samtools index -@ ${task.cpus} ${bam_id}.sorted.bam
    """
}

process FASTA_IDX {
    tag "$fasta"
    label "SAMTOOLS_CONTAINER"
    publishDir "${params.outdir}/fasta_idx", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta.baseName}.fa.fai", emit: fai

    script:
    """
    samtools faidx $fasta
    """
}

// PROCESS 1: Create Sequence Dictionary (Required for Picard)
process CREATE_DICT {
    tag "$fasta"
    label "PICARD_CONTAINER"

    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict

    script:
    """
    picard CreateSequenceDictionary \
        R=$fasta \
        O=${fasta.baseName}.dict
    """
}

// PROCESS 2: Convert BEDs to Interval Lists
process BEDS_TO_INTERVALS {
    tag "$probes_bed"
    tag "$targets_bed"
    label "PICARD_CONTAINER"

    input:
    path probes_bed
    path targets_bed
    path dict

    output:
    path "probes.interval_list", emit: bait_intervals // Renamed for clarity
    path "targets.interval_list", emit: target_intervals // Renamed for clarity

    script:
    """
    echo "Converting Probes BED to Interval List..."
    picard BedToIntervalList \
        I=$probes_bed \
        O=probes.interval_list \
        SD=$dict

    echo "Converting Targets BED to Interval List..."
    picard BedToIntervalList \
        I=$targets_bed \
        O=targets.interval_list \
        SD=$dict
    """
}

// PROCESS 3: Alignment Metrics (Chimeras, Unmapped %)
process PICARD_ALIGNMENT_METRICS {
    tag "$bam_id"
    label "PICARD_CONTAINER"

    publishDir "${params.outdir}/picard_alignment", mode: 'copy'

    input:
    tuple val(bam_id), path(bam)
    path ref

    output:
    path "*_alignment_metrics.txt", emit: metrics

    script:
    """
    picard CollectAlignmentSummaryMetrics \
        R=$ref \
        I=$bam \
        O=${bam_id}_alignment_metrics.txt
    """
}

// PROCESS 4: HS Metrics (On-Target Coverage, Fold 80 Penalty)
process PICARD_HS_METRICS {
    tag "$bam_id"
    label "PICARD_CONTAINER"

    publishDir "${params.outdir}/picard_hs", mode: 'copy'

    input:
    tuple val(bam_id), path(bam)
    path ref
    path idx
    path bait_intervals    // Now uses bait interval list
    path target_intervals  // Now uses target interval list

    output:
    path "*_hs_metrics.txt", emit: metrics

    script:
    """
    picard CollectHsMetrics \
        I=$bam \
        O=${bam_id}_hs_metrics.txt \
        R=$ref \
        BAIT_INTERVALS=$bait_intervals \
        TARGET_INTERVALS=$target_intervals
    """
}

// PROCESS 5: Estimate FPR via Ti/Tv (Bcftools)
process ESTIMATE_FPR_TITV {
    tag "$bam_id"
    publishDir "${params.outdir}/bcftools_stats", mode: 'copy'

    input:
    tuple val(bam_id), path(bam)
    path ref

    output:
    path "*_stats.txt", emit: stats

    script:
    """
    # 1. Mpileup, 2. Call variants, 3. Calculate Stats
    bcftools mpileup -Ou -f $ref $bam | \
    bcftools call -mv -Ob -o variants.bcf

    bcftools stats variants.bcf > ${bam_id}_stats.txt
    """
}

// PROCESS 6: MultiQC (Aggregation)
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path align_files
    path hs_files
    path bcf_stats

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

// WORKFLOW
workflow {
    // Channel setup
    bams_ch = Channel.fromPath(params.bams)
        .map { file -> tuple(file.simpleName, file) }
    bams_ch.view()
    ref_ch = Channel.value(file(params.ref))
    probes_bed_ch = Channel.value(file(params.probes_bed)) // New channel for probes
    targets_bed_ch = Channel.value(file(params.targets_bed)) // New channel for targets
    FASTA_IDX(ref_ch)
    
    // 1. Prepare Reference Dictionary
    CREATE_DICT(ref_ch)

    // 2. Prepare Interval Lists for Probes and Targets
    BEDS_TO_INTERVALS(probes_bed_ch, targets_bed_ch, CREATE_DICT.out.dict)
    // sort BAM
    SAMTOOLS_SORT(bams_ch)
    // 3. Run Picard Alignment Metrics (Chimeras, Unmapped)
    PICARD_ALIGNMENT_METRICS(SAMTOOLS_SORT.out.bam, ref_ch)

    // 4. Run Picard HS Metrics (Fold 80, On Target)
    //    Now passing both bait and target interval lists
    PICARD_HS_METRICS(
        bams_ch,
        ref_ch,
        FASTA_IDX.out.fai,
        BEDS_TO_INTERVALS.out.bait_intervals,    // Use bait intervals
        BEDS_TO_INTERVALS.out.target_intervals   // Use target intervals
    )

    // 5. Run Bcftools (FPR / TiTv)
    ESTIMATE_FPR_TITV(SAMTOOLS_SORT.out.bam, ref_ch)

    // 6. Aggregate results
    MULTIQC(
        PICARD_ALIGNMENT_METRICS.out.metrics.collect(),
        PICARD_HS_METRICS.out.metrics.collect(),
        ESTIMATE_FPR_TITV.out.stats.collect()
    )

}
