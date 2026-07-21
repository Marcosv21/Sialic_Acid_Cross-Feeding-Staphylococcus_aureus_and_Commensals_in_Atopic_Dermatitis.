#!/bin/bash
# Dependencies: GTDB-Tk
# Installation: conda install -c conda-forge -c bioconda gtdbtk
                

eval "$(conda shell.bash hook)"
conda activate gtdbtk

SOURCE_BINS_DIR="/temporario2/17404478/PRJEB59406/assay_1/gtdb_input_renamed/batch_1"
INPUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/gtdb_input_all_bins/batch_1"
OUTPUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/gtdb_taxonomy/batch_1"

# Set GTDB-Tk data path (adjust if necessary) 
# export GTDBTK_DATA_PATH="/temporario2/17404478/miniconda3/envs/gtdbtk/share/gtdbtk-2.7.2/db"

mkdir -p "$INPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Symlink all .fa files to the input directory
# -name: avoid recreating links if they already exist
# -exec ln -sf: create symbolic links, overwrite if they exist
# -sf: symbolic link, force overwrite if link exists
find "$SOURCE_BINS_DIR" -name "*.fa" -exec ln -sf {} "$INPUT_DIR" \;

# Validate file count
NUM_GENOMES=$(ls "$INPUT_DIR"/*.fa | wc -l)
echo "Genomes found: $NUM_GENOMES"

if [ "$NUM_GENOMES" -eq 0 ]; then 
    echo "ERROR: No .fa files found in $SOURCE_BINS_DIR"
    exit 1
fi

echo "Running GTDB-Tk classification..."

# Memory optimization: pplacer_cpus=1 and skip_ani_screen to reduce RAM usage
# Run the classification workflow
gtdbtk classify_wf \
    --genome_dir "$INPUT_DIR" \
    --extension 'fa' \
    --out_dir "$OUTPUT_DIR" \
    --cpus 16 \
    --pplacer_cpus 4 
    #--skip_ani_screen
# genome_dir: Directory with input genomes
# extension: File extension of genome files
# out_dir: Output directory
# cpus: Number of CPUs to use
# pplacer_cpus: Number of CPUs for pplacer (set to 1 for lower RAM usage)
# skip_ani_screen: Skip ANI screening to save memory, use with caution, may affect accuracy 
echo "------------------------------------------------"
echo "Done. Results: $OUTPUT_DIR"