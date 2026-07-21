#!/bin/bash

BINS_DIR="/temporario2/17404478/PRJEB59406/assay_1/magscot/batch_1/magscot_bins"

OUTPUT="/temporario2/17404478/PRJEB59406/assay_1/batch_1/mag_sample_mapping.tsv"

echo -e "MAG\tDominant_Sample\tNum_Contigs" > "$OUTPUT"

for bin in "$BINS_DIR"/*.fa; do

    MAG=$(basename "$bin")

    SAMPLE_INFO=$(grep "^>" "$bin" \
        | sed 's/^>//' \
        | cut -d'k' -f1 \
        | sort \
        | uniq -c \
        | sort -nr \
        | head -1)

    COUNT=$(echo "$SAMPLE_INFO" | awk '{print $1}')
    SAMPLE=$(echo "$SAMPLE_INFO" | awk '{print $2}')

    echo -e "${MAG}\t${SAMPLE}\t${COUNT}" >> "$OUTPUT"

done

echo "Done: $OUTPUT"