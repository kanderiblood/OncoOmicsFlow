#!/bin/bash
THREADS=16
WORKDIR=$(pwd)
OUTDIR="${WORKDIR}/quantms_out"

FASTA_DB="/HS_proteome.fasta"

QUANTMS_BIN="/quantms/bin"
FileConverter="${QUANTMS_BIN}/FileConverter"
CometAdapter="${QUANTMS_BIN}/CometAdapter"
PeptideIndexer="${QUANTMS_BIN}/PeptideIndexer"
PSMFeatureExtractor="${QUANTMS_BIN}/PSMFeatureExtractor"
PercolatorAdapter="${QUANTMS_BIN}/PercolatorAdapter"
FeatureFinderIdentification="${QUANTMS_BIN}/FeatureFinderIdentification"
ProteinQuantifier="${QUANTMS_BIN}/ProteinQuantifier"

PERCOLATOR_BIN="/percolator"
