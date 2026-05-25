#!/bin/bash

# Set variables
DEDUP_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/step3_dedup"
OUTPUT_DIR="${DEDUP_DIR}/bed_fastqc"
BEDTOOLS_PATH="/home/zhilin/anaconda3/bin/bedtools"
FASTQC_PATH="/home/zhilin/anaconda3/bin/fastqc"

# Create output directory (if it does not exist)
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/fastqc

# Find all deduplicated BAM files, convert to BED format, and run FASTQC analysis
for BAM_FILE in ${DEDUP_DIR}/dedup_*.bam; do
    SAMPLE_NAME=$(basename ${BAM_FILE} .bam)
    
    # Convert BAM to BED
    BED_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.bed"
    ${BEDTOOLS_PATH} bamtobed -i ${BAM_FILE} > ${BED_FILE}

    # Run FASTQC analysis
    ${FASTQC_PATH} -t 10 -o ${OUTPUT_DIR}/fastqc ${BAM_FILE}
done

