#!/bin/bash
# Dependencies: Bowtie2
# Install:
#   conda install -c bioconda bowtie2

# Activate the Conda environment
eval "$(conda shell.bash hook)"
conda activate Bowtie2

# Define paths
FASTQ_PATH="/temporario2/17404478/PRJEB59406/assay_1/cleaned_reads/batch_1"
INDEX_PATH="/temporario2/17404478/PRJEB59406/assay_1/indexed_contigs/batch_1"
OUTPUT_PATH="/temporario2/17404478/PRJEB59406/assay_1/contigs_reads_aligned/batch_1"

# Create the output directory if it does not exist
mkdir -p "$OUTPUT_PATH"

# Loop to process all *_R1.fastq files in the FASTQ_PATH folder
for FILE1 in "$FASTQ_PATH"/*_R1.fastq; do
  # Get the base sample name (without the _filtered_aligned_R1.fastq suffix)
  BASENAME=$(basename "$FILE1" _R1.fastq)

  # Define the corresponding file for the _R2.fastq pair
  FILE2="${FASTQ_PATH}/${BASENAME}_R2.fastq"
  INDEX_FILE="${INDEX_PATH}/${BASENAME}_indexed"

  # Check if the _R2.fastq file exists before running the alignment
  if [[ -f "$FILE2" ]]; then
    echo "Processing: $BASENAME"

    # Alignment with Bowtie2
    bowtie2 -x "$INDEX_FILE" \
      -1 "$FILE1" \
      -2 "$FILE2" \
      --threads 8 \
      -S "$OUTPUT_PATH/${BASENAME}_aligned.sam"
  else
    echo "Warning: Corresponding file $FILE2 not found. Skipping $BASENAME."
  fi

done