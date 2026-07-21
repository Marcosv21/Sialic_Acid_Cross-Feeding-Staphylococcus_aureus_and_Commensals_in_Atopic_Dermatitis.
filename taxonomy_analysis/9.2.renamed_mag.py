import os
import pandas as pd
import shutil

mapping_file = "/temporario2/17404478/PRJEB59406/assay_1/mag_sample_mapping.tsv"
source_dir = "/temporario2/17404478/PRJEB59406/assay_1/filtered_bins_high_quality/batch_1"
output_dir = "/temporario2/17404478/PRJEB59406/assay_1/gtdb_input_renamed/batch_1"

os.makedirs(output_dir, exist_ok=True)
df = pd.read_csv(mapping_file, sep="\t")
for _, row in df.iterrows():

    old_name = row["MAG"]
    sample = row["Dominant_Sample"]

    old_base = old_name.replace(".fa", "")

    mag_number = old_base.split("_")[-1]

    new_name = f"{sample}_MAG_{mag_number}.fa"

    src = os.path.join(source_dir, old_name)

    dst = os.path.join(output_dir, new_name)

    if os.path.exists(src):
        shutil.copy2(src, dst)
        print(f"{old_name} -> {new_name}")