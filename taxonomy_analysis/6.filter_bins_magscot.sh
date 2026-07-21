#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate seqkit

WORKDIR="/temporario2/17404478/PRJEB59406/assay_1/magscot/batch_1"
CONCAT_FASTA="$WORKDIR/all_contigs.fa"
MAGSCOT_OUT="$WORKDIR/MAGScoT_Final_Atopica.refined.contig_to_bin.out"
CONTIG_LIST_DIR="$WORKDIR/contig_lists"
BINS_DIR="$WORKDIR/magscot_bins"

THREADS=8

[[ ! -f "$CONCAT_FASTA" ]] && {
    echo "ERROR: FASTA not found"
    exit 1
}

[[ ! -f "$MAGSCOT_OUT" ]] && {
    echo "ERROR: MAGScoT output not found"
    exit 1
}

mkdir -p "$CONTIG_LIST_DIR"
mkdir -p "$BINS_DIR"

# Clean previous files
rm -f "$CONTIG_LIST_DIR"/*.txt
rm -f "$BINS_DIR"/*.fa

echo "Creating contig lists per bin..."

tail -n +2 "$MAGSCOT_OUT" | awk '
{
    bin=$1
    contig=$2

    # remove path from bin name
    gsub(/^.*\//, "", bin)

    print contig >> "'"$CONTIG_LIST_DIR"'/"bin".txt"
}
'

echo "Contig lists created."

echo "Extracting bins..."

for file in "$CONTIG_LIST_DIR"/*.txt; do

    bin=$(basename "$file" .txt)

    echo "Processing: $bin"

    seqkit grep \
        --threads "$THREADS" \
        -f "$file" \
        "$CONCAT_FASTA" \
        -o "$BINS_DIR/${bin}.fa"

done

TOTAL_BINS=$(ls "$BINS_DIR"/*.fa 2>/dev/null | wc -l)

echo ""
echo "=========================================="
echo "DONE"
echo "Generated bins: $TOTAL_BINS"
echo "Output directory:"
echo "$BINS_DIR"
echo "=========================================="