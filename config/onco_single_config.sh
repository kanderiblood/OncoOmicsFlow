#!/bin/bash
# Configuration for Onco panel pipeline

# Общее
THREADS=8
WORKDIR="${WORKDIR:-$(pwd)}"
BASE_DIR="${BASE_DIR:-/data/onco}"
REF_DIR="${REF_DIR:-$BASE_DIR/GRCh37}"
TOOLS_DIR="${TOOLS_DIR:-$BASE_DIR/tools}"
DB_DIR="${DB_DIR:-$BASE_DIR/db_37}"

# Референс
REF_GENOME="${REF_GENOME:-$REF_DIR/hs37d5.fa}"
ONCO_BED="${ONCO_BED:-$REF_DIR/nano.regions.final.nochr.bed}"
ONCO_REG="${ONCO_REG:-$REF_DIR/nano.regions.final.nochr.bed_1.interval_list}"

# Инструменты
FASTQC="${TOOLS_DIR}/fastqc"
FASTP="fastp"
BWA="bwa"
SAMTOOLS="samtools"
GATK="${TOOLS_DIR}/gatk-4.5.0.0/gatk"
PICARD_JAR="${TOOLS_DIR}/picard.jar"
MANTA="${TOOLS_DIR}/manta-1.6.0/bin"
VARSCAN_JAR="${TOOLS_DIR}/VarScan.v2.3.9.jar"
CNVKIT="cnvkit.py"
VEP_DIR="${TOOLS_DIR}/ensembl-vep-release-113/.vep"

# Conda environments
CONDA_ONCO="onco"
CONDA_BL="bl"
CONDA_CNVKIT="cnvkit_env"


# Панельные параметры
TUMOR_SITE=9
TUMOR_PURITY=0.9
TUMOR_PLOIDY=2.0
ASSAY="WGS"
