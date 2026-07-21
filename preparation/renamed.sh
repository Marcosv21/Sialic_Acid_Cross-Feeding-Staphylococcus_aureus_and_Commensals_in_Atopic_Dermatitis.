#!/bin/bash

INPUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/batch_1"
OUTPUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/renamed_reads/batch_1"
SAMPLE_MAP="/temporario2/17404478/PRJEB59406/sample_map.csv"

mkdir -p "$OUTPUT_DIR"

# In some cases, you need alter "IFS" for '\t', or delete "IFS=' '"
# In all cases, the sample_map.csv altered the "control" for "HC" 
tail -n +2 "$SAMPLE_MAP" | while IFS=',' read -r ERR SAMPLE REST; do
    R1="${INPUT_DIR}/${ERR}_1.fastq.gz"
    R2="${INPUT_DIR}/${ERR}_2.fastq.gz"
    OUT_R1="${OUTPUT_DIR}/${SAMPLE}_1.fastq.gz"
    OUT_R2="${OUTPUT_DIR}/${SAMPLE}_2.fastq.gz"

    if [[ ! -f "$R1" ]]; then
        echo "[SKIP] $SAMPLE: $R1"
        continue
    fi

    echo "[COPY] $ERR -> $SAMPLE"

    cp "$R1" "$OUT_R1"

    if [[ -f "$R2" ]]; then
        cp "$R2" "$OUT_R2"
    fi

done