#!/bin/bash
# Dependencies: MEGAHIT
# Install:
#   conda install -c bioconda megahit
eval "$(conda shell.bash hook)"

# 3. Activate environment
echo "activate megahit..."
conda activate megahit

FASTQ_PATH="/temporario2/17404478/PRJNA46333_2/assay/cleaned_reads" 
# MERGED_PATH="/temporario2/17404478/PRJNA46333_2/assay/merged_reads" 
OUTPUT_PATH="/temporario2/17404478/PRJNA46333_2/assay/megahit_assemblies"

mkdir -p "$OUTPUT_PATH"

#Loop to proccess all files in FASTQ_PATH
for FILE1 in "$FASTQ_PATH"/*_R1.fastq; do
  BASENAME=$(basename "$FILE1" _R1.fastq)          

  FILE2="${FASTQ_PATH}/${BASENAME}_R2.fastq"        
  # FILE3="${MERGED_PATH}/${BASENAME}.merged.fastq"   
  DIR_NAME="${OUTPUT_PATH}/${BASENAME}"
  
  # Check if assembly already exists
  if [ -f "$DIR_NAME/contigs.fa" ]; then
    echo "Assembly already exists for $BASENAME, skipping..."
    continue
  fi

  # Check is the incomplete paste exists, remove and reassemble
  if [ -d "$DIR_NAME" ]; then
    echo "Incomplete assembly found for $BASENAME, removing and reassembling..."
    rm -rf "$DIR_NAME"
  fi
  #Assembling with MEGAHIT
megahit \
    -1 "$FILE1" \
    -2 "$FILE2" \
    -o "$DIR_NAME" \
    -m 0.8 \
    -t 20 \
    --min-contig-len 500 \
    --k-min 27 \
    --k-max 127 \
    --k-step 10 \
    --merge-level 20,0.95 \
    --prune-level 2   
    # -r "$FILE3" \ -> Opcional, if you have merged reads, you can include them as unpaired reads in the assembly. However, if the merge percentage is low, it may not significantly improve the assembly and could increase computational time. You can test with and without this option to see if it makes a difference in your specific dataset.              
  # -1: File with forward reads (R1)
  # -2: File with reverse reads (R2)
  # -r: File with merged reads (unpaired) - optional
  # -o: Output directory for the assembly results
  # -m: Maximum memory usage (0.8 = 80% of available RAM)
  # -t: Number of threads to use
  # --min-contig-len: Minimum length of contigs to report (default is 200)
  # --k-min: Minimum k-mer size (default is 21)
  # --k-max: Maximum k-mer size (default is 141)
  # --k-step: Step size for k-mer sizes (default is 10)
  # --merge-level: Merge contigs with similar sequences (default is 20,0.95) (test 30,0.98)
  # --prune-level: Prune low-coverage tips from the assembly graph (test 3)

done
