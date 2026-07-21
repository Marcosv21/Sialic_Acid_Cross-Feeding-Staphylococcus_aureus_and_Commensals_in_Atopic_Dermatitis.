#!/bin/bash
# Dependências: Prodigal, Diamond
# Instalação: conda install prodigal diamond

# Inicializa o conda corretamente dentro do script
# source /home/marcos/miniconda3/etc/profile.d/conda.sh

# Agora você pode ativar o ambiente pelo nome
conda activate diamond
eval "$(conda shell.bash hook)"
conda activate diamond

MAGS_DIR="/temporario2/17404478/PRJEB59406/assay_1/gtdb_input_renamed/batch_1"

OUTPUT_DIR="/temporario2/17404478/PRJEB59406/assay_1/mag_annotation/batch_1"
PROTEINS_DIR="$OUTPUT_DIR/proteins"
DIAMOND_OUT="$OUTPUT_DIR/diamond_matches"

DB_PATH="/temporario2/17404478/PRJEB59406/code/data_base/pathway_db/db/sialic_acid_pathway.dmnd"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$PROTEINS_DIR"
mkdir -p "$DIAMOND_OUT"

for MAG in "$MAGS_DIR"/*.fa; do
    BIN_NAME=$(basename "$MAG" .fa)
    
    echo "Processing MAG: $BIN_NAME"
    prodigal -i "$MAG" \
             -a "$PROTEINS_DIR/${BIN_NAME}.faa" \
             -o "$PROTEINS_DIR/${BIN_NAME}.gff" \
             -p meta -q
             # -i: input fasta file
             # -a: output protein translations
             # -o: output GFF file
             # -p meta: metagenomic mode
             # -q: quiet mode (suppress output)

    diamond blastp -d "$DB_PATH" \
                   -q "$PROTEINS_DIR/${BIN_NAME}.faa" \
                   -o "$DIAMOND_OUT/${BIN_NAME}_hits_sial.tsv" \
                   --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                   --sensitive \
                   --quiet
                   # -d: Defining Diamond database path
                   # -q: Defining query file
                   # -o: Defining output file
                   # --outfmt 6: Define output format (tabular with specific fields)
                   # --mid-sensitive: Use mid-sensitive mode (balances speed and sensitivity)
                   # --quiet: Suppress standard output messages

    # Verification of hits
    if [ -s "$DIAMOND_OUT/${BIN_NAME}_hits_sial.tsv" ]; then
        echo " Found sialidase in $BIN_NAME!"
    fi
done
