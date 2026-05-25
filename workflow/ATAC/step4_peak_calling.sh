#!/bin/bash

# Set variables
BED_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/step3_dedup/bed_fastqc"
OUTPUT_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/step4_macs2"
GENOME_SIZE="1.35e8"  # Genome size for Arabidopsis thaliana
THREADS=10

# Create output directory (if it does not exist)
mkdir -p ${OUTPUT_DIR}

# Find all BED files and run MACS2 peak calling
for BED_FILE in ${BED_DIR}/*.bed; do
    SAMPLE_NAME=$(basename ${BED_FILE} .bed)
    macs2 callpeak -t ${BED_FILE} -f BEDPE -n ${SAMPLE_NAME} -g ${GENOME_SIZE} \
        --outdir ${OUTPUT_DIR} -q 0.05\
        --keep-dup all -B --SPMR --trackline
done

