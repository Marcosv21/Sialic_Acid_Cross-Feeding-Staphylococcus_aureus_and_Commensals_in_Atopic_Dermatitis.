#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate checkm

BINS_DIR="/temporario2/17404478/PRJEB59406/assay_1/magscot/batch_1/magscot_bins"
OUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/checkm_results/batch_1"

mkdir -p "$OUT_DIR"

count_bins=$(ls "$BINS_DIR"/*.fa 2>/dev/null | wc -l)

if [ "$count_bins" -eq 0 ]; then
    echo "ERROR: No bins found in $BINS_DIR"
    exit 1
fi

echo "Running CheckM on $count_bins MAGScoT-refined bins..."

checkm lineage_wf \
    -t 12 \
    --pplacer_threads 4 \
    -x fa \
    "$BINS_DIR" \
    "$OUT_DIR"

# Gera relatórios
if [ -f "$OUT_DIR/lineage.ms" ]; then
    echo "Generating QA reports..."

    checkm qa \
        "$OUT_DIR/lineage.ms" \
        "$OUT_DIR" \
        > "$OUT_DIR/quality_report.txt"

    checkm qa \
        "$OUT_DIR/lineage.ms" \
        "$OUT_DIR" \
        --tab_table \
        -f "$OUT_DIR/quality_table.tsv"

    echo "Done."
else
    echo "ERROR: lineage.ms not found — CheckM lineage_wf may have failed."
fi