#!/bin/bash

# Параметры потоков и директорий
THREADS=16
WORKDIR="${WORKDIR:-$(pwd)}"
OUTDIR="${WORKDIR}/results"

# Референсы
REF_DIR="${REF_DIR:-$WORKDIR/ref}"
REF_GENOME="${REF_GENOME:-$REF_DIR/GRCh37/hs37d5.fa}"
GTF_FILE="${GTF_FILE:-$REF_DIR/GRCh37/genes.gtf}"

# Панель для таргетных онко-панелей
ONCO_BED="${ONCO_BED:-$REF_DIR/panels/onco_panel.targets.bed}"
ONCO_PRO="${ONCO_PRO:-$REF_DIR/panels/onco_panel.probes.interval_list}"
ONCO_REG="${ONCO_REG:-$REF_DIR/panels/onco_panel.targets.interval_list}"

# Для CNV/TMB/MSI/Signatures
CNV_REF="${CNV_REF:-$REF_DIR/GRCh37/nano.flat_reference.cnn}"
MSI_DIR="${MSI_DIR:-$WORKDIR/msi}"
SIGNATURE_DIR="${SIGNATURE_DIR:-$WORKDIR/signatures}"

# Пути к инструментам
TOOLS_DIR="${TOOLS_DIR:-$WORKDIR/tools}"
FASTQC="${TOOLS_DIR}/fastqc"
FASTP="${TOOLS_DIR}/fastp"
BWA="${TOOLS_DIR}/bwa"
SAMTOOLS="${TOOLS_DIR}/samtools"
PICARD="${TOOLS_DIR}/picard.jar"
GATK="${TOOLS_DIR}/gatk-4.5.0.0/gatk"
MOSDEPTH="${TOOLS_DIR}/mosdepth"
CNVKIT="${TOOLS_DIR}/cnvkit.py"
MUTECT2="${GATK} Mutect2"
VEP_DIR="${TOOLS_DIR}/ensembl-vep-release-113"
