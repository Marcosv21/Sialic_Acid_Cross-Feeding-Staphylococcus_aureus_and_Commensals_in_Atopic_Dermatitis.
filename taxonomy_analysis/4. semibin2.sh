#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate semibin

# =========================
# CONFIG
# =========================
BASE_DIR="/temporario2/17404478/PRJEB59406/assay_1"
OUTPUT_DIR="$BASE_DIR/semibin2/batch_1"
BAM_DIR="$OUTPUT_DIR/realigned_bams"
TEMP_LINKS="$OUTPUT_DIR/temp_links"
FEATURE_DIR="$OUTPUT_DIR/features"
MODEL_DIR="$OUTPUT_DIR/models"
FINAL_DIR="$OUTPUT_DIR/final_bins"

MIN_READS=1000000
THREADS=12
MIN_CONTIG=1500

TOTAL_SAMPLES=0
VALID_SAMPLES=0
TRAINED_SAMPLES=0
BINNED_SAMPLES=0
Single_BINNED_SAMPLES=0
mkdir -p "$BAM_DIR" "$TEMP_LINKS" "$FEATURE_DIR" "$MODEL_DIR" "$FINAL_DIR"

# =========================
# STEP 1: CONCATENATE FASTA
# =========================
echo "Step 1: Preparing FASTA..."
rm -f "$TEMP_LINKS"/*.fa
for FOLDER in $BASE_DIR/megahit_assemblies/batch_1/*; do
    if [ -f "$FOLDER/final.contigs.fa" ]; then
        SAMPLE=$(basename "$FOLDER")
        seqkit seq -m $MIN_CONTIG "$FOLDER/final.contigs.fa" > "$TEMP_LINKS/${SAMPLE}.fa"

        N_CONTIGS=$(grep -c ">" "$TEMP_LINKS/${SAMPLE}.fa" 2>/dev/null || echo 0)
        if [[ "$N_CONTIGS" -eq 0 ]]; then
            echo "Warning: $SAMPLE has no contigs >= ${MIN_CONTIG}bp after filtering. Skipping."
            rm -f "$TEMP_LINKS/${SAMPLE}.fa"
            continue
        fi
        
        echo "  $SAMPLE: $N_CONTIGS contigs >= ${MIN_CONTIG}bp"
        ((TOTAL_SAMPLES++))
    fi
done

echo "Concatenating..."
SemiBin2 concatenate_fasta -i "$TEMP_LINKS"/*.fa -o "$OUTPUT_DIR/concatenated_dir"

FINAL_FASTA="$OUTPUT_DIR/concatenated.fa"
zcat "$OUTPUT_DIR/concatenated_dir"/*.fa.gz > "$FINAL_FASTA"

# STEP 2: INDEX
echo "Step 2: Indexing..."
bowtie2-build "$FINAL_FASTA" "$OUTPUT_DIR/idx"

# STEP 3: MAPPING
echo "Step 3: Mapping..."

VALID_BAMS=""
for FOLDER in $BASE_DIR/megahit_assemblies/batch_1/*; do
    SAMPLE=$(basename "$FOLDER")
    READ_DIR=$(dirname "$FOLDER")/../cleaned_reads

    R1="$READ_DIR/${SAMPLE}_R1.fastq"
    R2="$READ_DIR/${SAMPLE}_R2.fastq"

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "Mapping $SAMPLE..."

        BAM="$BAM_DIR/${SAMPLE}.bam"

        bowtie2 -p $THREADS \
            -x "$OUTPUT_DIR/idx" \
            -1 "$R1" -2 "$R2" 2> "$BAM_DIR/${SAMPLE}.log" | \
        samtools view -bS - | \
        samtools sort -o "$BAM"

        samtools index "$BAM"

        MAPPED=$(samtools view -c -F 4 "$BAM")

        if [[ "$MAPPED" -gt "$MIN_READS" ]]; then
            VALID_BAMS+="$BAM "
            ((VALID_SAMPLES++))
        else
            echo " $SAMPLE removed ($MAPPED reads)"
            rm -f "$BAM" "$BAM.bai"
        fi
    fi
done

if [ -z "$VALID_BAMS" ]; then
    echo "Error: no valid BAM files found"
    exit 1
fi

# STEP 4: FEATURES
echo "Step 4: Generating features..."

SemiBin2 generate_sequence_features_multi \
    --input-fasta "$FINAL_FASTA" \
    --input-bam $VALID_BAMS \
    --output "$FEATURE_DIR"

# STEP 5: TRAINING
echo "Step 5: Training models..."

for SAMPLE_PATH in "$FEATURE_DIR"/samples/*/; do

    SAMPLE=$(basename "$SAMPLE_PATH")

    DATA="$SAMPLE_PATH/data.csv"
    SPLIT="$SAMPLE_PATH/data_split.csv"

    if [[ -f "$DATA" && -f "$SPLIT" ]]; then

        N_LINES=$(wc -l < "$DATA")

        if [[ "$N_LINES" -gt 1 ]]; then
            echo "Training $SAMPLE ($N_LINES linhas)..."
            ((TRAINED_SAMPLES++))

            SemiBin2 train_self \
                --data "$DATA" \
                --data-split "$SPLIT" \
                --output "$MODEL_DIR/$SAMPLE" \
                --epochs 30

        else
            echo "Ignoring $SAMPLE (data.csv empty or insufficient)"
        fi

    else
        echo "Ignoring $SAMPLE (arquivos ausentes)"
    fi

done

# STEP 6: BINNING
echo "Step 6: Binning..."

for SAMPLE_PATH in "$FEATURE_DIR"/samples/*/; do
    SAMPLE=$(basename "$SAMPLE_PATH")

    INPUT_FASTA="$FEATURE_DIR/samples/${SAMPLE}.fa"
    MODEL="$MODEL_DIR/$SAMPLE/model.pt"
    DATA="$SAMPLE_PATH/data.csv"

    if [[ -f "$INPUT_FASTA" && -f "$MODEL" && -s "$DATA" ]]; then
        echo "Binning $SAMPLE..."

        SemiBin2 bin_short \
            --input-fasta "$INPUT_FASTA" \
            --model "$MODEL" \
            --data "$DATA" \
            --output "$FINAL_DIR/$SAMPLE"
        ((BINNED_SAMPLES++))
    else
        echo "Ignoring $SAMPLE in binning (incomplete data)"
        echo "Multi-sample failed → fallback to single-sample"
        mkdir -p "$OUTPUT_DIR/single.mode"
        SPECIFIC_BAM="$BAM_DIR/${SAMPLE}.bam"
        SPECIFIC_FASTA="$TEMP_LINKS/${SAMPLE}.fa"

        if [[ -f "$SPECIFIC_BAM" ]]; then
            SAMPLE=$(basename "$SPECIFIC_BAM" .bam)

            SemiBin2 single_easy_bin \
                --environment human_skin \
                -i "$SPECIFIC_FASTA" \
                -b "$SPECIFIC_BAM" \
                -o "$OUTPUT_DIR/single.mode/$SAMPLE" \
                > "$OUTPUT_DIR/single.mode/${SAMPLE}.log" 2>&1 || \
                echo "Failed: $SAMPLE"
            ((Single_BINNED_SAMPLES++))
        fi
    fi

done
# CLEANUP
rm -rf "$TEMP_LINKS"

echo "========================="
echo "Results:"
echo "Total samples: $TOTAL_SAMPLES"
echo "Samples with sufficient reads: $VALID_SAMPLES"
echo "Trained models: $TRAINED_SAMPLES"
echo "Samples with binning executed: $BINNED_SAMPLES"
echo "Samples with fallback to single-sample: $Single_BINNED_SAMPLES"
echo "========================="