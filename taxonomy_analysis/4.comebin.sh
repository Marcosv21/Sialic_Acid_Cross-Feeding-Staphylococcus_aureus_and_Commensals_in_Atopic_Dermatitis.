#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate comebin_env

MEGAHIT_DIR="/temporario2/17404478/PRJEB59406/assay_1/megahit_assemblies/batch_1"
BAM_DIR="/temporario2/17404478/PRJEB59406/assay_1/ordened_bams/batch_1"
OUTPUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/comebin_bins/batch_1"
THREADS=20
MIN_CONTIGS_FOR_COMEBIN=1000

mkdir -p "$OUTPUT_DIR"

TOTAL=0
DONE=0
SKIPPED=0
FAILED=0

for SAMPLE_DIR in "$MEGAHIT_DIR"/*/; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    CONTIGS_FILE="$SAMPLE_DIR/final.contigs.fa"
    BAM_FILE="$BAM_DIR/${SAMPLE_NAME}_aligned_sorted.bam"
    ((TOTAL++))
    echo "Sample: $SAMPLE_NAME"

    if [[ ! -f "$CONTIGS_FILE" ]]; then
        echo "[ERROR] Assembly not found:"
        echo "$CONTIGS_FILE"
        ((FAILED++))
        continue
    fi

    if [[ ! -f "$BAM_FILE" ]]; then
        echo "[ERROR] BAM not found:"
        echo "$BAM_FILE"
        ((FAILED++))
        continue
    fi

    TOTAL_CONTIGS=$(grep -c "^>" "$CONTIGS_FILE")

    echo "Total contigs: $TOTAL_CONTIGS"

    if [[ "$TOTAL_CONTIGS" -lt "$MIN_CONTIGS_FOR_COMEBIN" ]]; then
        echo "[SKIP] Too few contigs (< $MIN_CONTIGS_FOR_COMEBIN)"
        ((SKIPPED++))
        continue
    fi
    SAMPLE_OUTPUT="$OUTPUT_DIR/$SAMPLE_NAME"
    RESULT_BINS="$SAMPLE_OUTPUT/comebin_res/comebin_res_bins"
    mkdir -p "$SAMPLE_OUTPUT"

   if [[ -d "$RESULT_BINS" ]]; then
        N_EXISTING=$(find "$RESULT_BINS" \
            -name "*.fa" 2>/dev/null | wc -l)

        if [[ "$N_EXISTING" -gt 0 ]]; then
            echo "[SKIP] Already completed ($N_EXISTING bins)"
            ((DONE++))
            continue
        fi
    fi

    if [[ ! -f "${BAM_FILE}.bai" ]]; then
        echo "Creating BAM index..."
        samtools index \
            -@ "$THREADS" \
            "$BAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "[ERROR] Failed indexing BAM"
            ((FAILED++))
            continue
        fi
    fi
    BAM_LINK_DIR="$SAMPLE_OUTPUT/bam_input"

    rm -rf "$BAM_LINK_DIR"
    mkdir -p "$BAM_LINK_DIR"
    ln -sf "$BAM_FILE" \
        "$BAM_LINK_DIR/$(basename "$BAM_FILE")"

    ln -sf "${BAM_FILE}.bai" \
        "$BAM_LINK_DIR/$(basename "$BAM_FILE").bai"

    echo "Running COMEBin..."
    run_comebin.sh \
        -a "$CONTIGS_FILE" \
        -o "$SAMPLE_OUTPUT" \
        -p "$BAM_LINK_DIR" \
        -t "$THREADS" \
        -n 6

    STATUS=$?

    if [[ "$STATUS" -eq 0 ]] && [[ -d "$RESULT_BINS" ]]; then
        N_BINS=$(find "$RESULT_BINS" \
            -name "*.fa" | wc -l)
        echo "[SUCCESS] Generated $N_BINS bins"
        ((DONE++))

    else
        echo "[FAILED] COMEBin failed"
        ((FAILED++))

    fi
    rm -rf "$BAM_LINK_DIR"
done

echo "Total samples : $TOTAL"
echo "Completed     : $DONE"
echo "Skipped       : $SKIPPED"
echo "Failed        : $FAILED"
