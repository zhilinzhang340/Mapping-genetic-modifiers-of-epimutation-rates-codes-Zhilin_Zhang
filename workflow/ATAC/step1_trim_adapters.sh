#!/bin/bash

# Set variables
FASTQ_DIR="/mnt/ssd/RIL_meth/ATAC/raw_fastq"
OUTPUT_DIR="/mnt/ssd/RIL_meth/ATAC/pipeline_output/step1_trimmed"
THREADS=10  # Number of threads used by Trimmomatic
TRIMMOMATIC_DIR="/home/zhilin/documents/software/Trimmomatic-0.39"
TRIMMOMATIC_JAR="$TRIMMOMATIC_DIR/trimmomatic-0.39.jar"
ADAPTERS_FILE="$TRIMMOMATIC_DIR/adapters/TruSeq3-PE.fa"

# Create output directory (if it does not exist)
mkdir -p ${OUTPUT_DIR}

# Find all matching files and run adapter trimming
for SAMPLE_R1 in ${FASTQ_DIR}/*_1.fastq.gz; do
    SAMPLE_R2=${SAMPLE_R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE_NAME=$(basename ${SAMPLE_R1} _1.fastq.gz)

    # Trim adapter sequences
    java -jar ${TRIMMOMATIC_JAR} PE -threads ${THREADS} ${SAMPLE_R1} ${SAMPLE_R2} \
                    ${OUTPUT_DIR}/${SAMPLE_NAME}_R1_paired.fastq.gz ${OUTPUT_DIR}/${SAMPLE_NAME}_R1_unpaired.fastq.gz \
                    ${OUTPUT_DIR}/${SAMPLE_NAME}_R2_paired.fastq.gz ${OUTPUT_DIR}/${SAMPLE_NAME}_R2_unpaired.fastq.gz \
                    ILLUMINACLIP:${ADAPTERS_FILE}:1:30:9 \
                    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36
done

