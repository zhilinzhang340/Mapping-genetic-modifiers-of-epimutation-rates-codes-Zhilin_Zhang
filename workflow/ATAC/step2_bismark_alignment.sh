#!/bin/bash

# Set variables
GENOME_FILE="/home/zhilin/documents/new_ref/TAIR10_chr_all.fa"
OUTPUT_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/step2_alignment"
THREADS=10  # Number of threads for BWA and Samtools
SAMTOOLS_PATH="/home/zhilin/anaconda3/bin/samtools"

# Create output directory (if it does not exist)
mkdir -p ${OUTPUT_DIR}

# Find all matching files, and run alignment, sorting, indexing, and alignment statistics
for SAMPLE_R1 in /mnt/ssd/RIL_meth/ATAC/pipeline_output/step1_trimmed/*_R1_paired.fastq.gz; do
    SAMPLE_R2=${SAMPLE_R1/_R1_paired.fastq.gz/_R2_paired.fastq.gz}
    SAMPLE_NAME=$(basename ${SAMPLE_R1} _R1_paired.fastq.gz)

    # Perform alignment using BWA-MEM and directly convert and sort to BAM
    bwa mem -t ${THREADS} ${GENOME_FILE} ${SAMPLE_R1} ${SAMPLE_R2} | \
    ${SAMTOOLS_PATH} view -@ ${THREADS} -bS | \
    ${SAMTOOLS_PATH} sort -@ ${THREADS} -o ${OUTPUT_DIR}/sorted_${SAMPLE_NAME}.bam
    
    # Index the BAM file
    ${SAMTOOLS_PATH} index -@ ${THREADS} ${OUTPUT_DIR}/sorted_${SAMPLE_NAME}.bam

    # Generate alignment statistics
    ${SAMTOOLS_PATH} flagstat ${OUTPUT_DIR}/sorted_${SAMPLE_NAME}.bam > ${OUTPUT_DIR}/${SAMPLE_NAME}.stat
    cat ${OUTPUT_DIR}/${SAMPLE_NAME}.stat
done

