#!/bin/bash

# ================= 配置部分 =================

# 1. 输入与输出路径配置
INPUT_DIR="/mnt/ssd/RIL_meth/new_RNA-seq_OX_2026-05/raw_data/total_fastq"
OUTPUT_ROOT="/mnt/ssd/RIL_meth/new_RNA-seq_OX_2026-05/analysis_results"

# 2. 软件路径
TRIM_JAR="/home/zhilin/documents/software/Trimmomatic-0.38/trimmomatic-0.38.jar"
ADAPTERS="/home/zhilin/documents/software/Trimmomatic-0.38/adapters/TruSeq3-PE.fa"

# 3. 参考基因组路径
STAR_REF="/home/zhilin/documents/new_ref/STAR"
GTF_FILE="/home/zhilin/documents/ref/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.47.gtf"

# 4. 线程数
THREADS=12

# ================= 创建输出目录 =================
echo "正在创建输出目录: $OUTPUT_ROOT ..."

mkdir -p "${OUTPUT_ROOT}/00_merged_fastq"
mkdir -p "${OUTPUT_ROOT}/00_fastqc_raw"
mkdir -p "${OUTPUT_ROOT}/01_trimmed"
mkdir -p "${OUTPUT_ROOT}/02_aligned"
mkdir -p "${OUTPUT_ROOT}/03_counts"

# ================= 自动合并或链接 Fastq 文件 =================
echo "========================================================"
echo "开始检查并处理 Fastq 文件 (智能合并或软链接)..."
echo "读取目录: $INPUT_DIR"
echo "========================================================"

# 提取唯一的样本前缀 (适配 OX_Col_VIM... 格式，剔除 _23M..._L... 后缀)
SAMPLES=$(ls ${INPUT_DIR}/*_1.fq.gz 2>/dev/null | awk -F'/' '{print $NF}' | sed 's/_[^_]*_L[0-9]*_1\.fq\.gz//' | sort -u)

if [ -z "$SAMPLES" ]; then
    echo "错误: 在 $INPUT_DIR 中没有找到匹配的 *_1.fq.gz 文件，请检查路径。"
    exit 1
fi

for SAMPLE in $SAMPLES
do
    MERGED_R1="${OUTPUT_ROOT}/00_merged_fastq/${SAMPLE}_merged_1.fq.gz"
    MERGED_R2="${OUTPUT_ROOT}/00_merged_fastq/${SAMPLE}_merged_2.fq.gz"
    
    if [ ! -f "$MERGED_R1" ]; then
        # 将匹配到的文件存入 Bash 数组，方便统计数量
        R1_FILES=($(ls ${INPUT_DIR}/${SAMPLE}_*_1.fq.gz | sort))
        R2_FILES=($(ls ${INPUT_DIR}/${SAMPLE}_*_2.fq.gz | sort))
        FILE_COUNT=${#R1_FILES[@]}

        if [ "$FILE_COUNT" -gt 1 ]; then
            echo "样本 ${SAMPLE} 包含 $FILE_COUNT 个切片，正在执行 cat 合并..."
            cat "${R1_FILES[@]}" > "$MERGED_R1"
            cat "${R2_FILES[@]}" > "$MERGED_R2"
        elif [ "$FILE_COUNT" -eq 1 ]; then
            echo "样本 ${SAMPLE} 仅有 1 个切片，直接创建软链接 (节省空间与时间)..."
            # 使用绝对路径创建软链接，防止链接失效
            ln -s "${R1_FILES[0]}" "$MERGED_R1"
            ln -s "${R2_FILES[0]}" "$MERGED_R2"
        else
            echo "警告: 未找到样本 ${SAMPLE} 的文件！"
        fi
    else
        echo "样本 ${SAMPLE} 的处理文件已存在，跳过。"
    fi
done

echo "========================================================"
echo "原始数据准备完毕，开始执行核心分析流程..."
echo "========================================================"

# ================= 开始循环处理所有准备好的样本 =================

for r1_file in ${OUTPUT_ROOT}/00_merged_fastq/*_merged_1.fq.gz
do
    base_name=$(basename "$r1_file" "_1.fq.gz")
    r2_file="${OUTPUT_ROOT}/00_merged_fastq/${base_name}_2.fq.gz"
    
    # 获取纯净的样本 ID (例如 OX_Col_VIM2_Col_7-1)
    sample_id=${base_name%_merged}
    
    echo "########################################################"
    echo "正在处理样本: $sample_id"
    echo "使用线程数: $THREADS"
    echo "########################################################"

    # Step 0: FastQC
    echo "[Step 0] Running FastQC..."
    fastqc -o "${OUTPUT_ROOT}/00_fastqc_raw" -t $THREADS "$r1_file" "$r2_file"

    # Step 1: Trimmomatic
    echo "[Step 1] Running Trimmomatic..."
    java -jar $TRIM_JAR PE  -threads $THREADS -phred33 \
        "$r1_file" "$r2_file" \
        "${OUTPUT_ROOT}/01_trimmed/${sample_id}_1_paired.fq.gz"   "${OUTPUT_ROOT}/01_trimmed/${sample_id}_1_unpaired.fq.gz" \
        "${OUTPUT_ROOT}/01_trimmed/${sample_id}_2_paired.fq.gz"   "${OUTPUT_ROOT}/01_trimmed/${sample_id}_2_unpaired.fq.gz" \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Step 2: STAR Alignment
    echo "[Step 2] Running STAR Alignment..."
    STAR --genomeDir $STAR_REF \
         --runThreadN $THREADS \
         --readFilesIn "${OUTPUT_ROOT}/01_trimmed/${sample_id}_1_paired.fq.gz" "${OUTPUT_ROOT}/01_trimmed/${sample_id}_2_paired.fq.gz" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${OUTPUT_ROOT}/02_aligned/${sample_id}_" \
         --alignIntronMin 60 \
         --alignIntronMax 6000 \
         --outSAMtype BAM Unsorted

    # Step 3: Samtools Sorting
    echo "[Step 3] Sorting BAM file..."
    samtools sort -@ $THREADS "${OUTPUT_ROOT}/02_aligned/${sample_id}_Aligned.out.bam" \
        -o "${OUTPUT_ROOT}/02_aligned/${sample_id}_sorted.bam"
    
    samtools index "${OUTPUT_ROOT}/02_aligned/${sample_id}_sorted.bam"
    
    # Step 4: HTSeq-count
    echo "[Step 4] Running HTSeq-count..."
    htseq-count --format=bam \
                --mode=union \
                --stranded=reverse \
                --order=pos \
                --type=exon \
                "${OUTPUT_ROOT}/02_aligned/${sample_id}_sorted.bam" \
                "$GTF_FILE" > "${OUTPUT_ROOT}/03_counts/${sample_id}_counts.txt"

    echo ">>> 样本 $sample_id 处理完成。"
    echo ""

done

echo "所有流程结束！"
echo "请检查结果目录: $OUTPUT_ROOT"
