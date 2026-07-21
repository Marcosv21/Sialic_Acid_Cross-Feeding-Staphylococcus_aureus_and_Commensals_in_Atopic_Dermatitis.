#If you not ready use magscot, do you need to rename the contigs in the bins to match GTDB-Tk requirements?

#!/bin/bash

# =========================================================
# Rename contig headers + convert .fa -> .fna
# =========================================================

BINS_DIR="/temporario2/17404478/PRJNA46333_2/assay/filtered_bins_high_quality_semibin2"

cd $BINS_DIR || exit

echo "Processing MAGs in:"
echo $BINS_DIR
echo "-----------------------------------"

for f in *.fa
do

    # nome do MAG
    base=$(basename "$f" .fa)

    echo "Processing: $f"

    awk -v b="$base" '
    /^>/ {
        gsub(/^>/, "", $0)
        print ">" b "_" $0
        next
    }
    {print}
    ' "$f" > "${base}.tmp"

    # remove original
    rm "$f"

    # salva como .fna
    mv "${base}.tmp" "${base}.fna"

    echo "  -> Created ${base}.fna"

done

echo "-----------------------------------"
echo "DONE!"