#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate coverm

# DIRECTORIES

MAG_DIR="/temporario2/17404478/PRJEB59406/assay_1/gtdb_input_renamed/batch_1"
READS_DIR="/temporario2/17404478/PRJEB59406/assay_1/cleaned_reads/batch_1"
WORK_DIR="/temporario2/17404478/PRJEB59406/assay_1/mag_abundance/batch_1"

BAM_DIR="${WORK_DIR}/bam"
INDEX_DIR="${WORK_DIR}/bowtie2_index"
COVERM_DIR="${WORK_DIR}/coverm"

mkdir -p "$BAM_DIR"
mkdir -p "$INDEX_DIR"
mkdir -p "$COVERM_DIR"

# CONCATENATE ALL MAGS

ALL_MAGS="${INDEX_DIR}/all_mags.fa"
echo "Concatenating MAGs..."
cat ${MAG_DIR}/*.fa > "$ALL_MAGS"

# BUILD BOWTIE2 INDEX

echo "Building Bowtie2 index..."

bowtie2-build \
    "$ALL_MAGS" \
    "${INDEX_DIR}/mags_db"

# ALIGN EACH SAMPLE

echo "Starting alignment..."

for R1 in ${READS_DIR}/*_R1*.fastq
do

    SAMPLE=$(basename "$R1" | sed 's/_R1.*//')
    R2="${READS_DIR}/${SAMPLE}_R2.fastq"

    echo "Processing sample: $SAMPLE"

    # SAM output
    SAM_FILE="${BAM_DIR}/${SAMPLE}.sam"

    # BAM output
    BAM_FILE="${BAM_DIR}/${SAMPLE}.sorted.bam"

     if [[ -f "$BAM_FILE" ]]; then
        echo "  [SKIP] BAM already exists for $SAMPLE"
        continue
    fi

    # ALIGN

    bowtie2 \
        -x "${INDEX_DIR}/mags_db" \
        -1 "$R1" \
        -2 "$R2" \
        -p 16 \
        -S "$SAM_FILE"

    # SAM -> SORTED BAM

    samtools view -@ 16 -bS "$SAM_FILE" | \
    samtools sort -@ 16 -o "$BAM_FILE"

    # INDEX BAM
    samtools index "$BAM_FILE"

    # REMOVE HUGE SAM
    rm "$SAM_FILE"

done

# RUN COVERM

echo "Running CoverM..."

coverm genome \
    --bam-files ${BAM_DIR}/*.sorted.bam \
    --genome-fasta-directory "$MAG_DIR" \
    --genome-fasta-extension fa \
    --methods count \
    --min-covered-fraction 0 \
    --threads 10 \
    --output-file "${COVERM_DIR}/mag_abundance.tsv"

echo "DONE!"
echo "Results:"
echo "${COVERM_DIR}/mag_abundance.tsv"