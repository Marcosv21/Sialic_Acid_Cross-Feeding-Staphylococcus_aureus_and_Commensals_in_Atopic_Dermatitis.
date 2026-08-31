import pandas as pd
import glob

# ── CAMINHOS ─────────────────────────────────────────────────────────────────
GTDB_FILES = glob.glob("/home/marcos/PRJNA489681/table/gtdb_taxonomy/gtdbtk.*.summary.tsv")
ABUND_FILE = "/home/marcos/PRJNA489681/table/mag_abundance/coverm/mag_abundance_renamed.tsv"
OUT_FILE   = "/home/marcos/PRJNA489681/table/species_abundance_all.tsv"

# 1. Taxonomy (MAG -> Species) ────────────────────────────────────────────
df_tax = pd.concat(pd.read_csv(f, sep="\t", usecols=['user_genome', 'classification']) for f in GTDB_FILES)

mag_ids = df_tax['user_genome'].str.replace(".fa", "", regex=False)
species = df_tax['classification'].str.extract(r's__([^;]+)')[0].fillna('Unknown')

tax_dict = dict(zip(mag_ids, species))

# ── 2. SUBSTITUIR E SOMAR CONTAGENS ──────────────────────────────────────────
df_abund = pd.read_csv(ABUND_FILE, sep="\t")

# Mapping MAG IDs to Species and aggregating abundances
df_species = (
    df_abund.assign(Genome=df_abund['Genome'].map(tax_dict).fillna('Unknown'))
    .query("Genome != 'Unknown'")
    .groupby('Genome', as_index=False).sum()
    .rename(columns={'Genome': 'Species'})
)

# 3. Save ────────────────────────────────────────────────────────────────
df_species.to_csv(OUT_FILE, sep='\t', index=False)
print(f"table saved to: {OUT_FILE}")