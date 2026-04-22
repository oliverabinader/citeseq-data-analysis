#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 OUTPUT_BASE_DIR TRANSCRIPTOME FEATURE_REF LIBRARY_DIR"
    exit 1
fi

OUTPUT_BASE_DIR="$1"
TRANSCRIPTOME="$2"
FEATURE_REF="$3"
LIBRARY_DIR="$4"

for LIBRARY_FILE in "$LIBRARY_DIR"/*.csv; do
    
    SAMPLE_ID=$(basename "$LIBRARY_FILE" .csv)
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_ID}"

    mkdir -p "$OUTPUT_DIR"

    MRO_DISK_SPACE_CHECK=disable cellranger count \
        --id="$SAMPLE_ID" \
        --transcriptome="$TRANSCRIPTOME" \
        --create-bam=true \
        --output-dir="$OUTPUT_DIR" \
        --libraries="$LIBRARY_FILE" \
        --feature-ref="$FEATURE_REF"

    echo "Completed ${SAMPLE_ID}"

done
