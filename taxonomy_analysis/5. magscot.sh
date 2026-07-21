#!/bin/bash

#If you need to install MAGScoT, you can do it with:
# conda create -n magscott_env -c conda-forge -c bioconda r-base r-dplyr r-tidyr r-magrittr r-getopt perl -y
# conda activate magscott_env
# cd /temporario2/17404478/code/taxonomy_functions/
# git clone https://github.com/ikmb/MAGScoT.git
# cd MAGScoT
# chmod +x MAGScoT.R
# conda install -c conda-forge r-optparse -y
# Rscript MAGScoT.R --help

eval "$(conda shell.bash hook)"
conda activate magscott_env

# =========================
# PATHS
# =========================
BASE="/temporario2/17404478/PRJEB59406"
ASSEMBLY_DIR="$BASE/assay_1"
WORKDIR="$BASE/magscot/batch_1"
MAGSCOT_SCRIPT="/temporario2/17404478/code/taxonomy_functions/MAGScoT/MAGScoT.R"
HMM_SOURCE="/temporario2/17404478/code/taxonomy_functions/MAGScoT/hmm"

mkdir -p "$WORKDIR"
CONCAT_FASTA="$WORKDIR/all_contigs.fa"

# =========================
# STEP 0 — DEPENDÊNCIAS
# =========================
echo "Checking R packages..."
Rscript -e "libs<-c('digest','funr','optparse','dplyr','tidyr','readr'); \
    lapply(libs, function(x) if(!require(x,character.only=T)) \
    stop(paste('Missing package:',x)))"

# =========================
# STEP 1 — CONCATENAR CONTIGS
# =========================
echo "Step 1: Concatenating contigs..."
> "$CONCAT_FASTA"

for f in "$ASSEMBLY_DIR"/megahit_assemblies/batch_1/*/final.contigs.fa; do
    sample=$(basename "$(dirname "$f")")
    echo "  Processing $sample"
    awk -v s="$sample" '
        /^>/ { print ">" s "_" substr($0,2); next }
        { print }
    ' "$f" | seqkit seq -m 1500 >> "$CONCAT_FASTA"
done

N_CONTIGS=$(grep -c ">" "$CONCAT_FASTA")
echo "  Total contigs: $N_CONTIGS"

# =========================
# STEP 2 — HMM MARKERS
# =========================
echo "Step 2: Preparing HMM markers..."
hmmconvert "$HMM_SOURCE/gtdbtk_rel207_Pfam-A.hmm"  > "$WORKDIR/Pfam_clean.hmm"
hmmconvert "$HMM_SOURCE/gtdbtk_rel207_tigrfam.hmm" > "$WORKDIR/Tigrfam_clean.hmm"
cat "$WORKDIR/Pfam_clean.hmm" "$WORKDIR/Tigrfam_clean.hmm" > "$WORKDIR/bac120_markers.hmm"
hmmstat "$WORKDIR/bac120_markers.hmm" | tail -n 3

# =========================
# STEP 3 — PRODIGAL
# =========================
echo "Step 3: Running Prodigal..."
prodigal \
    -i "$CONCAT_FASTA" \
    -a "$WORKDIR/proteins.faa" \
    -p meta \
    > /dev/null 2>&1

echo "  Proteins predicted: $(grep -c ">" "$WORKDIR/proteins.faa")"

# =========================
# STEP 4 — HMMER
# =========================
echo "Step 4: Running hmmsearch..."
hmmsearch \
    --tblout "$WORKDIR/markers.tsv" \
    --cut_ga \
    --cpu 12 \
    "$WORKDIR/bac120_markers.hmm" \
    "$WORKDIR/proteins.faa"

awk '!/^#/ {print $1"\t"$3"\t"$13}' "$WORKDIR/markers.tsv" \
    > "$WORKDIR/magscot_hmm.tsv"

echo "  HMM hits: $(wc -l < "$WORKDIR/magscot_hmm.tsv")"

# =========================
# STEP 5 — BIN-CONTIG MAPPING
# =========================
echo "Step 5: Building bin-contig mapping..."
MAP_FILE="$WORKDIR/bin_mapping.tsv"
> "$MAP_FILE"

# --- MetaBAT2 ---
echo "  -> Collecting MetaBAT2 bins..."
METABAT_COUNT=0
for f in "$BASE"/assay_1/MetaBAT2_bins/batch_1/*/final.contigs.fa.metabat-bins*/*.fa; do
    [[ -f "$f" ]] || continue
    bin_name=$(basename "$f" .fa)
    sample=$(basename "$(dirname "$(dirname "$f")")")
    grep "^>" "$f" | sed 's/^>//' | \
    awk -v s="$sample" -v b="$bin_name" \
        '{print s"_"b"\t"s"_"$1"\tmetabat2"}' \
    >> "$MAP_FILE"
    ((METABAT_COUNT++))
done
echo "     MetaBAT2 bins: $METABAT_COUNT"

# --- SemiBin2 ---
echo "  -> Collecting SemiBin2 bins..."
SEMIBIN_COUNT=0
for d in "$BASE"/assay_1/semibin2/batch_1/final_bins/*/; do
    sample=$(basename "$d")
    [[ -f "$d/contig_bins.tsv" ]] || continue
    awk -v s="$sample" \
        '{print s"_sb2_"$2"\t"s"_"$1"\tsemibin2"}' \
        "$d/contig_bins.tsv" >> "$MAP_FILE"
    ((SEMIBIN_COUNT++))
done
echo "     SemiBin2 samples: $SEMIBIN_COUNT"

# --- COMEBin ---
echo "  -> Collecting COMEBin bins..."
COMEBIN_COUNT=0
for sample_dir in "$BASE"/assay_1/comebin_bins/batch_1/*/; do
    sample=$(basename "$sample_dir")
    BIN_DIR="$sample_dir/comebin_res/comebin_res_bins"

    [[ -d "$BIN_DIR" ]] || continue

    for f in "$BIN_DIR"/*.fa; do
        [[ -f "$f" ]] || continue
        bin_name=$(basename "$f" .fa)
        grep "^>" "$f" | sed 's/^>//' | \
        awk -v s="$sample" -v b="$bin_name" \
            '{print s"_cb_"b"\t"s"_"$1"\tcomebin"}' \
        >> "$MAP_FILE"
        ((COMEBIN_COUNT++))
    done
done
echo "     COMEBin bins: $COMEBIN_COUNT"

TOTAL_BINS=$((METABAT_COUNT + SEMIBIN_COUNT + COMEBIN_COUNT))
echo "  Total bins in mapping: $TOTAL_BINS"
echo "  Total lines in MAP_FILE: $(wc -l < "$MAP_FILE")"

# =========================
# STEP 6 — MAGSCOT
# =========================
echo "Step 6: Running MAGScoT..."

Rscript "$MAGSCOT_SCRIPT" \
    -i "$MAP_FILE" \
    --hmm "$WORKDIR/magscot_hmm.tsv" \
    -p bac120 \
    -o "$WORKDIR/MAGScoT_Final_Atopica" \
    -t 0.5 \
    --max_cont 0.15

echo "Done."
echo "Output: $WORKDIR/MAGScoT_Final_Atopica.*"