#!/bin/bash
# Dependencies: Samtools
# Install:
#   conda install -c bioconda samtools

#Activate conda environment
eval "$(conda shell.bash hook)"
conda activate samtools

# Path to the folder containing .sam files
SAM_DIR="/temporario2/17404478/PRJEB59406/assay_1/contigs_reads_aligned/batch_1"
# Path to save the sorted .bam files
OUTPUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/ordened_bams/batch_1"

mkdir -p $OUTPUT_DIR

for SAM_FILE in $SAM_DIR/*.sam; do
    # Extract the base filename (without extension)
    BASENAME=$(basename "$SAM_FILE" .sam)

    # Define names for intermediate and final files
    BAM_FILE="$OUTPUT_DIR/$BASENAME.bam"
    SORTED_BAM="$OUTPUT_DIR/${BASENAME}_sorted.bam"

    # Convert SAM to BAM
    echo "Converting $SAM_FILE to BAM..."
    samtools view -S -b "$SAM_FILE" > "$BAM_FILE"
    if [ $? -eq 0 ]; then
        echo "Removing $SAM_FILE..."
        rm -f "$SAM_FILE"
    else
        echo "Error converting $SAM_FILE. Keeping it for debugging."
        continue
    fi

     # Sort BAM
    echo "Sorting $BAM_FILE..."
    samtools sort "$BAM_FILE" -o "$SORTED_BAM"
    if [ $? -eq 0 ]; then
        echo "Removing intermediate BAM file $BAM_FILE..."

    else
        echo "Error sorting $BAM_FILE. Keeping the file for debugging."
    fi

    echo "Processing of $SAM_FILE completed!"
done

echo "All files have been successfully processed!"