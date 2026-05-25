#!/bin/bash

# Set variables
INPUT_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/step2_alignment"
OUTPUT_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/step3_dedup"
THREADS=10  # Number of threads used by Samtools
PICARD_PATH="/home/zhilin/documents/software/picard.jar"  # Path to Picard
SAMTOOLS_PATH="/home/zhilin/anaconda3/bin/samtools"

# Create output directory (if it does not exist)
mkdir -p ${OUTPUT_DIR}

# Find all matching BAM files and perform duplicate removal
for SAMPLE_BAM in ${INPUT_DIR}/sorted_*.bam; do
    SAMPLE_NAME=$(basename ${SAMPLE_BAM} .bam)

    # Remove duplicates
    java -jar ${PICARD_PATH} MarkDuplicates I=${SAMPLE_BAM} O=${OUTPUT_DIR}/dedup_${SAMPLE_NAME}.bam M=${OUTPUT_DIR}/${SAMPLE_NAME}_metrics.txt REMOVE_DUPLICATES=true

    # Index the deduplicated BAM file
    ${SAMTOOLS_PATH} index -@ ${THREADS} ${OUTPUT_DIR}/dedup_${SAMPLE_NAME}.bam

    # Generate alignment statistics for the deduplicated BAM
    ${SAMTOOLS_PATH} flagstat ${OUTPUT_DIR}/dedup_${SAMPLE_NAME}.bam > ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup.flagstat
    ${SAMTOOLS_PATH} idxstats ${OUTPUT_DIR}/dedup_${SAMPLE_NAME}.bam > ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup.idxstats
    cat ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup.flagstat
    cat ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup.idxstats
done

