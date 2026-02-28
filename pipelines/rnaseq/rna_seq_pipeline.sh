#!/bin/bash
# Скрипт для анализа данных RNA-seq


CONFIG_FILE="$(dirname "$0")/../../config/rna_config.sh"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Config file not found: $CONFIG_FILE"
    exit 1
fi
source "$CONFIG_FILE"


mkdir -p "${OUTDIR}/01_fastqc_raw" \
         "${OUTDIR}/02_trimmed" \
         "${OUTDIR}/03_align" \
         "${OUTDIR}/04_quant" \
         "${OUTDIR}/05_qc_summary"


for sample_line in "${SAMPLES[@]}"; do
    read -r SAMPLEID R1 R2 <<< "$sample_line"
    echo "Processing sample: $SAMPLEID"

    "$FASTQC" -t "$THREADS" -o "${OUTDIR}/01_fastqc_raw" "$R1" "$R2"

    TRIM_R1="${OUTDIR}/02_trimmed/${SAMPLEID}_R1_trimmed.fq.gz"
    TRIM_R2="${OUTDIR}/02_trimmed/${SAMPLEID}_R2_trimmed.fq.gz"

    "$FASTP" -i "$R1" -I "$R2" \
        -o "$TRIM_R1" -O "$TRIM_R2" \
        --detect_adapter_for_pe --thread "$THREADS" \
        --html "${OUTDIR}/02_trimmed/${SAMPLEID}_fastp.html" \
        --json "${OUTDIR}/02_trimmed/${SAMPLEID}_fastp.json"

    "$FASTQC" -t "$THREADS" -o "${OUTDIR}/02_trimmed" "$TRIM_R1" "$TRIM_R2"
done


"$MULTIQC" "${OUTDIR}/01_fastqc_raw" -o "${OUTDIR}/05_qc_summary/multiqc_raw"
"$MULTIQC" "${OUTDIR}/02_trimmed" -o "${OUTDIR}/05_qc_summary/multiqc_trimmed"


for sample_line in "${SAMPLES[@]}"; do
    read -r SAMPLEID R1 R2 <<< "$sample_line"
    SAMPLE_DIR="${OUTDIR}/03_align/${SAMPLEID}"
    mkdir -p "$SAMPLE_DIR"

    BAM_FILE="${SAMPLE_DIR}/${SAMPLEID}_Aligned.sortedByCoord.out.bam"

    "$STAR_BIN" \
        --runThreadN "$THREADS" \
        --genomeDir "$REF_GENOME_STAR_INDEX" \
        --readFilesIn "${OUTDIR}/02_trimmed/${SAMPLEID}_R1_trimmed.fq.gz" \
                       "${OUTDIR}/02_trimmed/${SAMPLEID}_R2_trimmed.fq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_DIR}/${SAMPLEID}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMattributes NH HI AS nM XS

    "$SAMTOOLS" index "$BAM_FILE"
    "$SAMTOOLS" flagstat "$BAM_FILE" > "${SAMPLE_DIR}/${SAMPLEID}_flagstat.txt"
done


MERGED_RSEM="${OUTDIR}/05_qc_summary/rsem_tpm_merged.txt"
MERGED_FC="${OUTDIR}/05_qc_summary/featureCounts_merged.txt"

> "$MERGED_RSEM"
> "$MERGED_FC"

HEADER_RSEM_DONE=false
HEADER_FC_DONE=false

for sample_line in "${SAMPLES[@]}"; do
    read -r SAMPLEID R1 R2 <<< "$sample_line"
    BAM_FILE="${OUTDIR}/03_align/${SAMPLEID}/${SAMPLEID}_Aligned.sortedByCoord.out.bam"
    QUANT_DIR="${OUTDIR}/04_quant/${SAMPLEID}"
    mkdir -p "$QUANT_DIR"

    "$RSEM_BIN" \
        --paired-end \
        --bam \
        --estimate-rspd \
        --calc-pme \
        --calc-ci \
        --num-threads "$THREADS" \
        "$BAM_FILE" "$REF_GENOME_RSEM_INDEX" \
        "${QUANT_DIR}/${SAMPLEID}"

    RSEM_TPM="${QUANT_DIR}/${SAMPLEID}.genes.results"
    if [ "$HEADER_RSEM_DONE" = false ]; then
        awk 'NR==1{print "GeneID",$6} NR>1{print $1,$6}' "$RSEM_TPM" > "$MERGED_RSEM"
        HEADER_RSEM_DONE=true
    else
        awk 'NR>1{print $6}' "$RSEM_TPM" | paste "$MERGED_RSEM" - > "${MERGED_RSEM}.tmp"
        mv "${MERGED_RSEM}.tmp" "$MERGED_RSEM"
    fi

    FC_FILE="${QUANT_DIR}/${SAMPLEID}_featureCounts.txt"
    "$FEATURECOUNTS_BIN" -T "$THREADS" -p -B -C -a "$GTF_FILE" -o "$FC_FILE" "$BAM_FILE"

    if [ "$HEADER_FC_DONE" = false ]; then
        awk 'NR==1{print $1,$7} NR>1{print $1,$7}' "$FC_FILE" > "$MERGED_FC"
        HEADER_FC_DONE=true
    else
        awk 'NR>1{print $7}' "$FC_FILE" | paste "$MERGED_FC" - > "${MERGED_FC}.tmp"
        mv "${MERGED_FC}.tmp" "$MERGED_FC"
    fi
done


Z_SCORE_FILE="${OUTDIR}/05_qc_summary/rsem_tpm_zscore.txt"
Rscript -e "
t <- read.table('$MERGED_RSEM', header=TRUE, row.names=1)
z <- t
for (i in 1:nrow(t)) {
    z[i,] <- (t[i,] - mean(t[i,])) / sd(t[i,])
}
write.table(z, file='$Z_SCORE_FILE', sep='\t', quote=FALSE, col.names=NA)
"
