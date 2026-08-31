import pandas as pd
import glob
import os

# ── CAMINHOS ─────────────────────────────────────────────────────────────────
GTDB_FILES = glob.glob("/home/marcos/PRJNA489681/table/gtdb_taxonomy/gtdbtk.*.summary.tsv")
DIAMOND_FILES = glob.glob("/home/marcos/PRJNA489681/table/mag_annotation/diamond_matches/*_hits_sial.tsv")
OUT_FILE = "/home/marcos/PRJNA489681/table/species_diamond_matches.tsv"

GENES = ['nanH', 'nanE', 'nanK', 'nanA', 'nanT']

# ── 1. Taxonomy (MAG -> Species) ────────────────────────────────────────────
df_tax = pd.concat(pd.read_csv(f, sep="\t", usecols=['user_genome', 'classification']) for f in GTDB_FILES)
tax_dict = dict(zip(
    df_tax['user_genome'].str.replace(".fa", "", regex=False),
    df_tax['classification'].str.extract(r's__([^;]+)')[0].fillna('Unknown')
))

# ── 2. Read Diamond Files and Check Genes ────────────────────────────────────────────
mag_data = []
for f in DIAMOND_FILES:
    mag_id = os.path.basename(f).replace("_hits_sial.tsv", "")
    if os.path.getsize(f) > 0:
        # Pega a coluna de matches e junta tudo num texto só para busca rápida
        hits = pd.read_csv(f, sep="\t", header=None)[1].str.cat(sep=' ').lower()
        mag_data.append({
            'mag_id': mag_id,
            **{f'has_{g}': (g.lower() in hits) for g in GENES}
        })

df_diamond = pd.DataFrame(mag_data)

# ── 3. Substitute MAG IDs with Species and Group ──────────────────────────────────────
df_diamond['Species'] = df_diamond['mag_id'].map(tax_dict).fillna('Unknown')
df_diamond = df_diamond.query("Species != 'Unknown'").drop(columns=['mag_id'])

df_species = df_diamond.groupby('Species', as_index=False).any()

for g in GENES:
    df_species[f'has_{g}'] = df_species[f'has_{g}'].map({True: 'YES', False: 'NO'})

# ── 4. Save ────────────────────────────────────────────────────────────────
df_species.to_csv(OUT_FILE, sep='\t', index=False)
print(f"table saved to: {OUT_FILE}")