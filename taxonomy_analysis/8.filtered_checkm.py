import os
import pandas as pd
import shutil
import glob


CHECKM_DIR = "/temporario2/17404478/PRJEB59406/assay_1/checkm_results/batch_1"
BINS_DIR = "/temporario2/17404478/PRJEB59406/assay_1/magscot/batch_1/magscot_bins"
OUTPUT_DIR = "/temporario2/17404478/PRJEB59406/assay_1/filtered_bins_high_quality/batch_1"

# Quality Thresholds (MIMAG Medium)
#If you want High Quality bins, set MIN_COMP = 90.0 and MAX_CONT = 5.0. High quality bins are suitable for most downstream analyses.
#For Medium Quality bins, set MIN_COMP = 50.0 and MAX_CONT = 10.0. Medium quality bins are generally acceptable for most downstream analyses.
#For Low Quality bins, set MIN_COMP = 0.0 and MAX_CONT = 25.0. Low quality bins are generally not recommended for downstream analyses.
#Article reference: Bowers et al., 2017. Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. Nat Biotechnol 35, 725–731. https://doi.org/10.1038/nbt.3893
MIN_COMP = 50.0
MAX_CONT = 10.0

# ---------------------
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"--- Filtering: Completeness >= {MIN_COMP}%, Contamination <= {MAX_CONT}% ---")
tables = glob.glob(f"{CHECKM_DIR}/quality_table.tsv")
total_copied = 0
for table in tables:
    try:
        # Load CheckM table
        df = pd.read_csv(table, sep="\t")
        # Filter good bins
        good_bins = df[
            (df['Completeness'] >= MIN_COMP) &
            (df['Contamination'] <= MAX_CONT)
        ]
        if good_bins.empty:
            continue
        print(f"Processing {len(good_bins)} bins to copy.")
        # Copy bins
        for bin_id in good_bins['Bin Id']:
            bin_file = os.path.join(BINS_DIR, f"{bin_id}.fa")
            if os.path.isfile(bin_file):
                dst = os.path.join(OUTPUT_DIR, f"{bin_id}.fa")
                shutil.copy2(bin_file, dst)
                total_copied += 1
                print(f"  Copied: {bin_id}.fa")
            else:
                print(f"  Warning: File not found for {bin_id}")
    except Exception as e:
        print(f"Error: {e}")
print("------------------------------------------------")
print(f"DONE! {total_copied} bins copied to:")
print(OUTPUT_DIR)