#!/bin/bash
#Dependencies: fastp, fastqc
#Install:
#  conda install -c bioconda fastp fastqc

eval "$(conda shell.bash hook)"
conda activate fastp

# Define paths
INPUT_DIR="/temporario2/17404478/PRJEB63144/assay/renamed_reads"
OUTPUT_DIR="/temporario2/17404478/PRJEB63144/assay/fastp_filtered"
REPORT_DIR="/temporario2/17404478/PRJEB63144/assay/fastp_quality"
FASTQC_DIR="/temporario2/17404478/PRJEB63144/assay/fastqc_results"

# Create output directories if they do not exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORT_DIR"
mkdir -p "$FASTQC_DIR"

FAILED_LOG="${OUTPUT_DIR}/fastp_failed.txt"
> "$FAILED_LOG"

for FILE1 in "$INPUT_DIR"/*_1.fastq.gz; do
    BASENAME=$(basename "$FILE1" "_1.fastq.gz")
    FILE2="${INPUT_DIR}/${BASENAME}_2.fastq.gz"

    OUTPUT1="${OUTPUT_DIR}/${BASENAME}_filtered_1.fastq.gz"
    OUTPUT2="${OUTPUT_DIR}/${BASENAME}_filtered_2.fastq.gz"
    HTML_REPORT="${REPORT_DIR}/${BASENAME}_fastp.html"
    JSON_REPORT="${REPORT_DIR}/${BASENAME}_fastp.json"

    if [[ ! -f "$FILE2" ]]; then
        echo "Warning: Corresponding file $FILE2 not found. Skipping $BASENAME." | tee -a "$FAILED_LOG"
        continue
    fi

    if ! gzip -t "$FILE1" 2>/dev/null || ! gzip -t "$FILE2" 2>/dev/null; then
        echo "Warning: Corrupted gzip for $BASENAME. Skipping." | tee -a "$FAILED_LOG"
        continue
    fi

    echo "-------------------------------------------------------"
    echo "Processing $BASENAME..."
    
    # 1. Executa o fastp
    fastp -i "$FILE1" -I "$FILE2" \
          -o "$OUTPUT1" -O "$OUTPUT2" \
          --detect_adapter_for_pe \
          -h "$HTML_REPORT" \
          -j "$JSON_REPORT"

    if [ $? -ne 0 ]; then
        echo "Error processing fastp for $BASENAME!" | tee -a "$FAILED_LOG"
        continue
    fi

    echo "Running FastQC for $BASENAME..."
    fastqc -o "$FASTQC_DIR" "$OUTPUT1" "$OUTPUT2"

    if [ $? -ne 0 ]; then
        echo "Error running FastQC for $BASENAME!" | tee -a "$FAILED_LOG"
    fi

done

echo "-------------------------------------------------------"
echo "Processing complete!"
echo "fastp reports: $REPORT_DIR"
echo "FastQC reports: $FASTQC_DIR"

FAIL_COUNT=$(wc -l < "$FAILED_LOG")
if [ "$FAIL_COUNT" -gt 0 ]; then
    echo "Some files were skipped. Check $FAILED_LOG"
fi