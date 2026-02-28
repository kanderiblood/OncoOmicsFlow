#!/bin/bash
# Скрипт для анализа одного образца с QuantMS (без расчета статистики)

set -euo pipefail

if [[ -z "$RAW_FILE" ]]; then
    echo "Укажите входной файл (.raw или .mzML)"
    echo "Пример: bash quantms_single_sample.sh /path/to/sample.raw"
    exit 1
fi

if [[ ! -f "$FASTA_DB" ]]; then
    echo "FASTA база не найдена: $FASTA_DB"
    exit 1
fi

SAMPLE=$(basename "$RAW_FILE" | sed 's/\.[^.]*$//')
OUTDIR="${OUTDIR}/${SAMPLE}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

if [[ "$RAW_FILE" == *.raw ]]; then
    echo "Converting RAW to mzML..."
    "$FileConverter" -in "$RAW_FILE" -out "${SAMPLE}.mzML"
    IN_MZML="${SAMPLE}.mzML"
else
    IN_MZML="$RAW_FILE"
fi

"$CometAdapter" \
    -in "$IN_MZML" \
    -database "$FASTA_DB" \
    -out "${SAMPLE}_comet.idXML" \
    -enzyme1 Trypsin \
    -fragment_mass_tolerance 0.02 \
    -precursor_mass_tolerance 10.0 \
    -fixed_modifications "Carbamidomethyl (C)" \
    -variable_modifications "Oxidation (M)" \
    -max_variable_mods_per_peptide 3 \
    -max_precursor_charge 7 \
    -threads "$THREADS"

"$PeptideIndexer" \
    -in "${SAMPLE}_comet.idXML" \
    -fasta "$FASTA_DB" \
    -out "${SAMPLE}_indexed.idXML" \
    -decoy_string DECOY_ \
    -threads "$THREADS"


"$PSMFeatureExtractor" \
    -in "${SAMPLE}_indexed.idXML" \
    -out "${SAMPLE}_feat.idXML"


echo "Running Percolator for FDR..."
"$PercolatorAdapter" \
    -in "${SAMPLE}_feat.idXML" \
    -out "${SAMPLE}_perc.idXML" \
    -percolator_executable "$PERCOLATOR_BIN/percolator" \
    -enzyme trypsin \
    -decoy_pattern DECOY_ \
    -score_type pep \
    -threads "$THREADS" \
    -protein_level_fdr_cutoff 0.05


"$FeatureFinderIdentification" \
    -in "$IN_MZML" \
    -id "${SAMPLE}_perc.idXML" \
    -out "${SAMPLE}.featureXML"


echo "Quantifying proteins..."
"$ProteinQuantifier" \
    -in "${SAMPLE}.featureXML" \
    -ids "${SAMPLE}_perc.idXML" \
    -protein_groups_out "${SAMPLE}_proteins.tsv" \
    -peptide_out "${SAMPLE}_peptides.tsv" \
    -top 3 -average

echo "QuantMS processing for sample $SAMPLE finished."
