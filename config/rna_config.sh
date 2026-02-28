#!/bin/bash
THREADS=16
WORKDIR=$(pwd)
OUTDIR="${WORKDIR}/rna_seq_out"

FASTQC="${WORKDIR}/tools/fastqc/fastqc"
FASTP="${WORKDIR}/tools/fastp/fastp"
MULTIQC="${WORKDIR}/tools/multiqc/multiqc"

STAR_BIN="${WORKDIR}/tools/STAR/STAR"
SAMTOOLS="${WORKDIR}/tools/samtools/samtools"
RSEM_BIN="${WORKDIR}/tools/rsem/rsem-calculate-expression"
FEATURECOUNTS_BIN="${WORKDIR}/tools/subread/featureCounts"

REF_GENOME_STAR_INDEX="${WORKDIR}/db/hg38/star_index"
REF_GENOME_RSEM_INDEX="${WORKDIR}/db/hg38/rsem_index"
GTF_FILE="${WORKDIR}/db/hg38/annotations/Homo_sapiens.GRCh38.109.gtf"

SAMPLES=(
    "Sample1 ${WORKDIR}/reads/sample1_R1.fastq.gz ${WORKDIR}/reads/sample1_R2.fastq.gz"
    "Sample2 ${WORKDIR}/reads/sample2_R1.fastq.gz ${WORKDIR}/reads/sample2_R2.fastq.gz"
)


#MULTIQC_CONFIG="${WORKDIR}/tools/multiqc/multiqc_config.yaml"
