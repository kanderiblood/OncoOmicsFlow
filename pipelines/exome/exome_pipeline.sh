#!/bin/bash

TUMOR_SAMPLE=$1
NORMAL_SAMPLE=$2

if [[ -z "$TUMOR_SAMPLE" || -z "$NORMAL_SAMPLE" ]]; then
    echo "Usage: $0 <TUMOR_SAMPLE> <NORMAL_SAMPLE>"
    exit 1
fi


source config/exome_config.sh

mkdir -p "$WORKDIR"/{fastq,pre_align_qc,trim,align,post_align_qc,calling/{somatic/{mutect,vardict,varscan,strelka,manta,cnv},germline},vep}


SAMPLES=$(ls "$FASTQ_DIR" | grep "_R1_" | sed 's/_R1_.*.fastq.gz//' | sort | uniq)

for SAMPLE in $SAMPLES; do
    echo ">>> FASTQ QC & trim: $SAMPLE"

    R1="$FASTQ_DIR/${SAMPLE}_R1.fastq.gz"
    R2="$FASTQ_DIR/${SAMPLE}_R2.fastq.gz"

    fastqc -t $THREADS -o "$WORKDIR/pre_align_qc" "$R1" "$R2"

    fastp \
        -i "$R1" -I "$R2" \
        -o "$WORKDIR/trim/${SAMPLE}_R1_trim.fq.gz" \
        -O "$WORKDIR/trim/${SAMPLE}_R2_trim.fq.gz" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --n_base_limit 5 \
        --thread $THREADS \
        --html "$WORKDIR/trim/${SAMPLE}_fastp.html" \
        --json "$WORKDIR/trim/${SAMPLE}_fastp.json"

    fastqc -t $THREADS -o "$WORKDIR/trim" "$WORKDIR/trim/${SAMPLE}_R1_trim.fq.gz" "$WORKDIR/trim/${SAMPLE}_R2_trim.fq.gz"
done

multiqc "$WORKDIR/trim" -o "$WORKDIR/trim/multiqc"


for SAMPLE in $SAMPLES; do
    echo ">>> Alignment: $SAMPLE"

    R1="$WORKDIR/trim/${SAMPLE}_R1_trim.fq.gz"
    R2="$WORKDIR/trim/${SAMPLE}_R2_trim.fq.gz"

    BAM_ALIGN="$WORKDIR/align/${SAMPLE}.bam"
    BAM_SORT="$WORKDIR/align/${SAMPLE}_sorted.bam"
    BAM_FINAL="$WORKDIR/post_align_qc/${SAMPLE}_final.bam"
    METRICS="$WORKDIR/post_align_qc/${SAMPLE}_metrics.txt"
    RECAL_TABLE="$WORKDIR/align/${SAMPLE}_recal.table"

    bwa mem -t $THREADS -Y -R "@RG\tID:$SAMPLE\tPL:ILLUMINA\tPU:$SAMPLE\tSM:$SAMPLE" \
        "$REF_GENOME" "$R1" "$R2" > "$WORKDIR/align/${SAMPLE}.sam"

    samtools view -@ $THREADS -bhS "$WORKDIR/align/${SAMPLE}.sam" > "$BAM_ALIGN"

    java -jar "$PICARD" SortSam I="$BAM_ALIGN" O="$BAM_SORT" SO=coordinate
    java -jar "$PICARD" MarkDuplicates I="$BAM_SORT" O="$BAM_FINAL" M="$METRICS" REMOVE_DUPLICATES=false

    samtools index "$BAM_FINAL"

    "$GATK" BaseRecalibrator -R "$REF_GENOME" -I "$BAM_FINAL" \
        --known-sites "$WORKDIR/db/dbsnp_138.b37.vcf.gz" \
        --known-sites "$WORKDIR/db/1000G_phase1.indels.b37.vcf.gz" \
        --known-sites "$WORKDIR/db/Mills_and_1000G_gold_standard.indels.b37.vcf.gz" \
        -O "$RECAL_TABLE" -L "$ONCO_REG"

    "$GATK" ApplyBQSR -R "$REF_GENOME" -I "$BAM_FINAL" --bqsr-recal-file "$RECAL_TABLE" -O "$WORKDIR/post_align_qc/${SAMPLE}.final.bqsr.bam"
    samtools index "$WORKDIR/post_align_qc/${SAMPLE}.final.bqsr.bam"
done

BAMS=$(ls "$WORKDIR/post_align_qc/*.final.bqsr.bam")
for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" .final.bqsr.bam)
    echo ">>> Post-alignment QC: $SAMPLE"

    java -jar "$PICARD" CollectInsertSizeMetrics \
        I="$BAM" \
        O="$WORKDIR/post_align_qc/${SAMPLE}_insert_size_metrics.txt" \
        H="$WORKDIR/post_align_qc/${SAMPLE}_insert_size_metrics.pdf"

    java -jar "$PICARD" CollectAlignmentSummaryMetrics \
        I="$BAM" \
        O="$WORKDIR/post_align_qc/${SAMPLE}_alignment_metrics.txt" \
        R="$REF_GENOME"

    "$GATK" CollectSequencingArtifactMetrics \
        -I "$BAM" \
        -O "$WORKDIR/post_align_qc/${SAMPLE}_artifact_metrics.txt" \
        -R "$REF_GENOME"

    java -jar "$PICARD" CollectHsMetrics \
        I="$BAM" \
        O="$WORKDIR/post_align_qc/${SAMPLE}_hs_metrics.txt" \
        R="$REF_GENOME" \
        BI="$ONCO_REG" \
        TI="$ONCO_PRO"

    mosdepth -t $THREADS -b "$ANOT_BED" "$SAMPLE" "$BAM"
done

multiqc "$WORKDIR/post_align_qc" -o "$WORKDIR/post_align_qc/multiqc"


# --- Mutect2 ---
MUTECT_RAW="$WORKDIR/calling/somatic/mutect/${TUMOR_SAMPLE}.raw.vcf.gz"
MUTECT_FILTERED="$WORKDIR/calling/somatic/mutect/${TUMOR_SAMPLE}.filtered.vcf.gz"

"$GATK" Mutect2 \
    -R "$REF_GENOME" \
    -I "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" -tumor "$TUMOR_SAMPLE" \
    -I "$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam" -normal "$NORMAL_SAMPLE" \
    -L "$ANOT_BED" \
    --germline-resource "$WORKDIR/db/af-only-gnomad.raw.sites.grch37.vcf.gz" \
    --f1r2-tar-gz "$WORKDIR/calling/somatic/mutect/f1r2.tar.gz" \
    -O "$MUTECT_RAW"

"$GATK" GetPileupSummaries \
    -I "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" \
    -V "$WORKDIR/db/gnomad.vcf.gz" -L "$ANOT_BED" \
    -O "$WORKDIR/calling/somatic/mutect/pileups.table"

"$GATK" CalculateContamination \
    -I "$WORKDIR/calling/somatic/mutect/pileups.table" \
    -O "$WORKDIR/calling/somatic/mutect/contamination.table"

"$GATK" LearnReadOrientationModel \
    -I "$WORKDIR/calling/somatic/mutect/f1r2.tar.gz" \
    -O "$WORKDIR/calling/somatic/mutect/read-orientation-model.tar.gz"

"$GATK" FilterMutectCalls \
    -V "$MUTECT_RAW" \
    --contamination-table "$WORKDIR/calling/somatic/mutect/contamination.table" \
    --ob-priors "$WORKDIR/calling/somatic/mutect/read-orientation-model.tar.gz" \
    -O "$MUTECT_FILTERED"

bcftools filter -i 'DP>=20 && INFO/AF>=0.05 && QUAL>=30 && FILTER="PASS"' \
    "$MUTECT_FILTERED" -Oz -o "$WORKDIR/calling/somatic/mutect/mutect.final.vcf.gz"

# --- VarScan2 ---
cd "$WORKDIR/calling/somatic/varscan"
samtools mpileup -l "$ANOT_BED" --no-BAQ -f "$REF_GENOME" --min-BQ 20 \
    "$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam" \
    "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" \
    | java -jar "$VARSCAN" somatic varscan --mpileup 1 --output-vcf

java -jar "$VARSCAN" processSomatic varscan.snp.vcf varscan.snp
java -jar "$VARSCAN" processSomatic varscan.indel.vcf varscan.indel
bgzip *.vcf && tabix *.vcf.gz

bcftools concat -a varscan.snp.Somatic.vcf.gz varscan.indel.Somatic.vcf.gz \
    | bcftools filter -i 'DP>=20 && FORMAT/AF>=0.05 && QUAL>=30' -Oz \
    -o "$WORKDIR/calling/somatic/varscan/varscan.final.vcf.gz"

# --- VarDict ---
cd "$WORKDIR/calling/somatic/vardict"
vardict-java -G "$REF_GENOME" -f 0.01 \
    -N "$TUMOR_SAMPLE" -c 1 -S 2 -E 3 -g 4 -th $THREADS \
    -b "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam|$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam" \
    "$ANOT_BED" | testsomatic.R | var2vcf_paired.pl \
    -N "$TUMOR_SAMPLE|$NORMAL_SAMPLE" -S -M -f 0.01 \
    > vardict.raw.vcf

bcftools filter -i 'DP>=20 && QUAL>=30 && FORMAT/AF>=0.05' \
    vardict.raw.vcf -Oz -o "$WORKDIR/calling/somatic/vardict/vardict.final.vcf.gz"

# --- Strelka2 ---
cd "$WORKDIR/calling/somatic/strelka"
configureStrelkaSomaticWorkflow.py \
    --normalBam "$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam" \
    --tumorBam "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" \
    --referenceFasta "$REF_GENOME" \
    --targeted \
    --runDir .

python2 runWorkflow.py -m local -j $THREADS

bcftools concat \
    results/variants/somatic.snvs.vcf.gz \
    results/variants/somatic.indels.vcf.gz \
    | bcftools filter -i 'FORMAT/DP>=20 && FORMAT/AF>=0.05' -Oz \
    -o "$WORKDIR/calling/somatic/strelka/strelka.final.vcf.gz"

# --- Manta SV ---
"$MANTA_DIR/configManta.py" \
    --normalBam "$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam" \
    --tumorBam "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" \
    --exome \
    --referenceFasta "$REF_GENOME" \
    --runDir "$WORKDIR/calling/somatic/manta/"

"$WORKDIR/calling/somatic/manta/runWorkflow.py" -m local -j $THREADS -g 60

cd "$WORKDIR/calling/somatic/cnv"
conda activate cnvkit

cnvkit.py autobin -f "$REF_GENOME" \
    -t "$ONCO_PRO" \
    --target-output-bed target.bed \
    --antitarget-output-bed antitarget.bed \
    "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" \
    "$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam"

cnvkit.py batch \
    "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" \
    --normal "$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam" \
    --fasta "$REF_GENOME" \
    --output-dir cnvkit_output/ \
    --scatter --diagram \
    -p $THREADS \
    --drop-low-coverage -t target.bed

TMB_FILE="$WORKDIR/somatic_metrics/${TUMOR_SAMPLE}_tmb.txt"

mkdir -p "$WORKDIR/somatic_metrics"

# фильтруем только exonic nonsynonymous
bcftools view -i 'INFO/EXONIC=="Y" && INFO/EFFECT~"missense_variant|nonsense_variant|frameshift_variant|splice_acceptor_variant|splice_donor_variant"' \
    "$WORKDIR/calling/somatic/mutect/mutect.final.vcf.gz" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
    | wc -l > "$TMB_FILE"


echo "scale=2; $(cat $TMB_FILE)/30" | bc > "$WORKDIR/somatic_metrics/${TUMOR_SAMPLE}_tmb_perMb.txt"

MSI_OUT="$WORKDIR/somatic_metrics/msi"

mkdir -p "$MSI_OUT"

msisensor scan -d "$REF_GENOME" -o "$MSI_OUT/microsatellites.list"

msisensor msi \
    -d "$MSI_OUT/microsatellites.list" \
    -n "$WORKDIR/post_align_qc/${NORMAL_SAMPLE}.final.bqsr.bam" \
    -t "$WORKDIR/post_align_qc/${TUMOR_SAMPLE}.final.bqsr.bam" \
    -o "$MSI_OUT/${TUMOR_SAMPLE}_msi"

python3 <<EOF
import pandas as pd
from SigProfilerExtractor import sigpro as sig

# input VCF (Mutect2 final)
input_vcf = "$WORKDIR/calling/somatic/mutect/mutect.final.vcf.gz"
output_dir = "$WORKDIR/somatic_metrics/signatures"

sig.sigProfilerExtractor(input_type="vcf", 
                         input_data=input_vcf, 
                         output=output_dir, 
                         reference_genome="GRCh37", 
                         exome=True, 
                         minimum_signatures=1, 
                         maximum_signatures=5)
EOF
echo ">>> Exome tumor-normal pipeline finished"
