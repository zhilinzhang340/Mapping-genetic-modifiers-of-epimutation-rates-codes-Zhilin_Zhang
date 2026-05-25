#!/bin/bash

# Set variables
FASTQ_DIR="/mnt/ssd/RIL_meth/ATAC/raw_fastq"
OUTPUT_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/fastqc"
THREADS=10 # Number of threads for each FastQC process

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Define quality control function
qc() {
    SAMPLE_R1=$1
    SAMPLE_R2=${SAMPLE_R1/_1.fastq.gz/_2.fastq.gz}

    # Run FastQC
    fastqc -t ${THREADS} ${SAMPLE_R1} ${SAMPLE_R2} -o ${OUTPUT_DIR}
}

export -f qc

# Find all matching files and perform quality control
find ${FASTQ_DIR} -name "*_1.fastq.gz" | while read SAMPLE_R1; do
    qc ${SAMPLE_R1}
done




