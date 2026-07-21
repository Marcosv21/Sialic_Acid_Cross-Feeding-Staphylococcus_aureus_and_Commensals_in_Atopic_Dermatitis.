# If you need process in different queues, use the command for concatenation different files in different queues, and then merge the results in a final step.
import os
import glob
import pandas as pd
from functools import reduce

base_dir = "/home/marcos/PRJEB59406/table/mag_abundance"

files = sorted(
    glob.glob(os.path.join(base_dir, "batch_*", "coverm", "mag_abundance.tsv"))
)

if len(files) == 0:
    raise FileNotFoundError("No mag_abundance.tsv files found.")

tables = []

for file in files:

    print(f"Reading {file}")

    df = pd.read_csv(file, sep="\t")

    # Rename H_ -> HC_
    new_columns = []

    for col in df.columns:

        if col.startswith("H_"):
            col = col.replace("H_", "HC_", 1)

        new_columns.append(col)

    df.columns = new_columns

    tables.append(df)

# Merge all tables by Genome
merged = reduce(
    lambda left, right: pd.merge(
        left,
        right,
        on="Genome",
        how="outer"
    ),
    tables
)

merged = merged.fillna(0)

# Convert counts to integers
for col in merged.columns[1:]:
    merged[col] = merged[col].astype(int)

output = os.path.join(base_dir, "mag_abundance_all.tsv")

merged.to_csv(
    output,
    sep="\t",
    index=False
)

print("\nDone!")
print(f"Output: {output}")
print(f"MAGs: {merged.shape[0]}")
print(f"Samples: {merged.shape[1]-1}")