#!/bin/bash
# Dependencies: Bowtie2, samtools
# Install:
#   conda install -c bioconda bowtie2 samtools
eval "$(conda shell.bash hook)"
conda activate bowtie2

FASTQ_PATH="/temporario2/17404478/PRJNA46333_2/assay/fastp_filtered/batch_1"
GENOME_INDEX_PATH="/temporario2/17404478/code/data_base/human_genome/GRCh38_index"
OUTPUT_PATH="/temporario2/17404478/PRJNA46333_2/assay/bowtie2_aligned/batch_1"

mkdir -p "$OUTPUT_PATH"

for FILE1 in "$FASTQ_PATH"/*_filtered_1.fastq.gz; do
  BASENAME=$(basename "$FILE1" "_filtered_1.fastq.gz")
  FILE2="${FASTQ_PATH}/${BASENAME}_filtered_2.fastq.gz"
  FILE_MERGED="${FASTQ_PATH}/${BASENAME}_merged.fastq.gz"

  echo "Processing: $BASENAME"

  # Paired-end — Filtered by Bowtie2 directly, then sorted and converted to BAM
  bowtie2 -x "$GENOME_INDEX_PATH" \
    -1 "$FILE1" \
    -2 "$FILE2" \
    --threads 8 \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    2> "$OUTPUT_PATH/${BASENAME}_paired.log" \
  | samtools view -@ 8 -b -f 12 -F 256 \
  | samtools sort -n -@ 8 \
    -o "$OUTPUT_PATH/${BASENAME}_filtered.bam"
# -x: path to the Bowtie2 index (without the .bt2 extension)
# -1: forward reads
# -2: reverse reads
# --threads: number of threads to use
# --very-sensitive: preset for high sensitivity
# --no-mixed: only report paired alignments (no singletons)
# --no-discordant: only report concordant alignments (no discordant pairs)
# 2>: redirect Bowtie2's stderr to a log file for troubleshooting

  echo "Done: $BASENAME"
done