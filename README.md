# OncoOmicsFlow 
Reproducible multi-omics workflow for clinical interpretation of WES, WGS, RNA-seq, targeted panels, and proteomics data. Modular bash pipelines, future migration to Nextflow + Docker planned.

# Features

-Proteomics: QuantMS, Comet, Percolator, ProteinQuantifier

-RNA-seq: FastQC, fastp, STAR, RSEM, featureCounts, MultiQC

-Exome / Targeted / WGS: BWA, Samtools, Picard, GATK (Mutect2, BaseRecalibrator), VarScan, Manta, CNVkit, Ensembl VEP, PCGR

Reproducibility: Modular bash scripts, multi-threaded, designed for cluster execution

Future-proof: Migration to Nextflow pipelines and Docker containers planned

# Example for proteomics
bash quantms_single_sample.sh /path/to/sample.raw

# Example for RNA-seq
bash rna_seq_pipeline.sh /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz

# Example for WES/Targeted tumor-only Requirements
bash onco_pipeline_tumor.sh SAMPLE_ID /path/to/R1.fastq.gz /path/to/R2.fastq.gz

# Requirements
-Linux environment with bash

-Installed dependencies as listed in each pipeline

-Access to reference genomes, GTF/GFF files, and database resources
