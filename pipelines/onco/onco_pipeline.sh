#!/bin/bash
# Скрипт только для единичного опухолевого образца

SAMPLE_ID="$1"
FASTQ_R1="$2"
FASTQ_R2="$3"

if [[ -z "$SAMPLE_ID" || -z "$FASTQ_R1" || -z "$FASTQ_R2" ]]; then
    echo "Usage: bash $0 SAMPLE_ID FASTQ_R1 FASTQ_R2"
    exit 1
fi

# Загружаем конфиг
CONFIG_FILE="$(dirname "$0")/../../config/onco_config.sh"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Config file not found: $CONFIG_FILE"
    exit 1
fi
source "$CONFIG_FILE"

START_TIME=$(date +%s)

cd "$WORKDIR"

# Создание каталогов
mkdir -p pre_align_qc trim post_align_qc align \
         calling/Mutect2 calling/varscan calling/cv \
         calling/cnv vep pcgr_outputs

# Активируем conda окружение 
source "$BASE_DIR/miniconda3/bin/activate"
conda activate "$CONDA_ONCO"

# QC и тримминг
"$FASTQC" -t "$THREADS" -o pre_align_qc "$FASTQ_R1" "$FASTQ_R2"

"$FASTP" -i "$FASTQ_R1" -I "$FASTQ_R2" \
         -o trim/${SAMPLE_ID}_R1_tr.fq \
         -O trim/${SAMPLE_ID}_R2_tr.fq \
         --detect_adapter_for_pe

"$FASTQC" -t "$THREADS" -o trim trim/${SAMPLE_ID}_R1_tr.fq trim/${SAMPLE_ID}_R2_tr.fq

# Выравнивание
"$BWA" mem -t "$THREADS" -Y \
    -R "@RG\tID:$SAMPLE_ID\tPL:ILLUMINA\tSM:$SAMPLE_ID" \
    "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" > align/${SAMPLE_ID}.sam

"$SAMTOOLS" view -@ "$THREADS" -bhS align/${SAMPLE_ID}.sam > align/${SAMPLE_ID}.bam
"$SAMTOOLS" index align/${SAMPLE_ID}.bam
"$SAMTOOLS" flagstat align/${SAMPLE_ID}.bam

# Picard: сортировка и дубликаты
java -jar "$PICARD_JAR" SortSam \
    I=align/${SAMPLE_ID}.bam \
    O=align/${SAMPLE_ID}_namesorted.bam \
    SO=queryname

java -jar "$PICARD_JAR" MarkDuplicates \
    I=align/${SAMPLE_ID}_namesorted.bam \
    O=post_align_qc/${SAMPLE_ID}_markdup.bam \
    ASSUME_SORT_ORDER=queryname \
    METRICS_FILE=post_align_qc/${SAMPLE_ID}_markdup_metrics.txt \
    VALIDATION_STRINGENCY=LENIENT

java -jar "$PICARD_JAR" SortSam \
    I=post_align_qc/${SAMPLE_ID}_markdup.bam \
    O=post_align_qc/${SAMPLE_ID}_sorted_markdup.bam \
    SO=coordinate

"$SAMTOOLS" index post_align_qc/${SAMPLE_ID}_sorted_markdup.bam

# GATK BaseRecalibrator и ApplyBQSR
"$GATK" BaseRecalibrator \
    -R "$REF_GENOME" \
    -I post_align_qc/${SAMPLE_ID}_sorted_markdup.bam \
    -O post_align_qc/${SAMPLE_ID}.recal.table \
    --known-sites "$DB_DIR/dbsnp_138.b37.vcf.gz" \
    --known-sites "$DB_DIR/1000G_phase1.indels.b37.vcf.gz" \
    --known-sites "$DB_DIR/Mills_and_1000G_gold_standard.indels.b37.vcf.gz" \
    -L "$ONCO_BED"

"$GATK" ApplyBQSR \
    -R "$REF_GENOME" \
    -I post_align_qc/${SAMPLE_ID}_sorted_markdup.bam \
    -O post_align_qc/${SAMPLE_ID}_bqsr.bam \
    --bqsr-recal-file post_align_qc/${SAMPLE_ID}.recal.table

# Picard Metrics
java -jar "$PICARD_JAR" CollectInsertSizeMetrics \
    I=post_align_qc/${SAMPLE_ID}_bqsr.bam \
    O=post_align_qc/${SAMPLE_ID}_insert_size_metrics.txt \
    H=post_align_qc/${SAMPLE_ID}_insert_size_metrics.pdf

java -jar "$PICARD_JAR" CollectAlignmentSummaryMetrics \
    I=post_align_qc/${SAMPLE_ID}_bqsr.bam \
    O=post_align_qc/${SAMPLE_ID}_alignment_metrics.txt \
    R="$REF_GENOME"

"$GATK" CollectSequencingArtifactMetrics \
    -I post_align_qc/${SAMPLE_ID}_bqsr.bam \
    -O post_align_qc/${SAMPLE_ID}_artifact_metrics.txt \
    -R "$REF_GENOME"

java -jar "$PICARD_JAR" CollectHsMetrics \
    I=post_align_qc/${SAMPLE_ID}_bqsr.bam \
    O=post_align_qc/${SAMPLE_ID}_hs_metrics.txt \
    R="$REF_GENOME" \
    BI="$ONCO_REG" \
    TI="$ONCO_PRO"

# Mosdepth
cd post_align_qc
mosdepth -t "$THREADS" -b "$ONCO_BED" ${SAMPLE_ID}_bqsr ${SAMPLE_ID}_bqsr.bam
cd ..

# FastQC BAM
"$FASTQC" -t "$THREADS" -o post_align_qc post_align_qc/${SAMPLE_ID}_bqsr.bam
multiqc post_align_qc

# Variant calling Mutect2
"$GATK" Mutect2 \
    -R "$REF_GENOME" \
    -I post_align_qc/${SAMPLE_ID}_bqsr.bam \
    --tumor-sample "$SAMPLE_ID" \
    -L "$ONCO_BED" \
    -O calling/Mutect2/raw.vcf \
    --panel-of-normals \
    --germline-resource "$DB_DIR/gnomad.vcf.gz"

"$GATK" FilterMutectCalls \
    -R "$REF_GENOME" \
    -V calling/Mutect2/raw.vcf \
    -O calling/Mutect2/filtered.vcf

bcftools filter -i 'FILTER="PASS"' calling/Mutect2/filtered.vcf > calling/Mutect2/pass.vcf

# VarScan
"$SAMTOOLS" mpileup -B -q 1 -f "$REF_GENOME" \
    -l "$ONCO_BED" post_align_qc/${SAMPLE_ID}_bqsr.bam > calling/varscan/tumor.mpileup

java -jar "$VARSCAN_JAR" mpileup2snp calling/varscan/tumor.mpileup \
    --output-vcf 1 --variants --min-var-freq 0.01 --min-coverage 10 > calling/varscan/snp.vcf

java -jar "$VARSCAN_JAR" mpileup2indel calling/varscan/tumor.mpileup \
    --output-vcf 1 --variants --min-var-freq 0.01 --min-coverage 10 > calling/varscan/indel.vcf

# Фильтрация VarScan
bcftools filter -i 'QUAL>=20 && INFO/DP>=20 && INFO/PVAL<=0.05 && INFO/AF>=0.02' calling/varscan/snp.vcf \
    > calling/varscan/snp.filtered.vcf

bcftools filter -i 'QUAL>=30 && INFO/DP>=30 && INFO/PVAL<=0.01 && INFO/AF>=0.05' calling/varscan/indel.vcf \
    > calling/varscan/indel.filtered.vcf

bcftools concat -a -O z -o calling/varscan/wgs.vcf.gz calling/varscan/snp.filtered.vcf calling/varscan/indel.filtered.vcf
bcftools merge -O z calling/varscan/wgs.vcf.gz calling/Mutect2/pass.vcf -o calling/merge_variants.vcf.gz

END_TIME=$(date +%s)
echo "Finished in $((END_TIME - START_TIME)) seconds"
