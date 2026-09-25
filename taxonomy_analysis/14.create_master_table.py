import pandas as pd
import os
import glob

# ── 1. PATHS ─────────────────────────────────────────────────────────────────
GTDB_DIR = "/home/marcos/PRJNA489681/table/gtdb_taxonomy"
DIAMOND_DIR = "/home/marcos/PRJNA489681/table/mag_annotation (2)/diamond_matches"
DESEQ2_FILE = "/home/marcos/PRJNA489681/table/deseq2_species_all_results.tsv"

OUTPUT_FILE = "/home/marcos/PRJNA489681/table/master_table_mag_sialidase_gtdb_deseq2.tsv"

GENES = ['has_nanH', 'has_nanE', 'has_nanK', 'has_nanA', 'has_nanT']

# ── FUNCTIONS ────────────────────────────────────────────────────────────────
def parse_taxonomy(tax_string):
    """Splits GTDB taxonomy into Phylum, Genus and Species"""
    if pd.isna(tax_string):
        return pd.Series(['Unknown', 'Unknown', 'Unknown'])
    
    parts = tax_string.split(';')
    p = next((x.split('__')[1] for x in parts if 'p__' in x), 'Unknown')
    g = next((x.split('__')[1] for x in parts if 'g__' in x), 'Unknown')
    s = next((x.split('__')[1] for x in parts if 's__' in x), 'Unknown')
    
    return pd.Series([p, g, s])

# ── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    
    # 1. LOAD TAXONOMY 
    tax_files = glob.glob(os.path.join(GTDB_DIR, "gtdbtk.*.summary.tsv"))
    if not tax_files:
        print("ERROR: No GTDB-Tk files found.")
        return
        
    df_tax = pd.concat([pd.read_csv(f, sep="\t", usecols=['user_genome', 'classification']) for f in tax_files])
    df_tax['mag_id'] = df_tax['user_genome'].str.replace(".fa", "", regex=False).str.replace(".fasta", "", regex=False)
    df_tax[['Phylum', 'Genus', 'Species']] = df_tax['classification'].apply(parse_taxonomy)
    
    df_tax = df_tax[['mag_id', 'Phylum', 'Genus', 'Species']]
    
    print(f"Taxonomy loaded: {len(df_tax)} MAGs")

    # 2. LOAD AND AGGREGATE DIAMOND FUNCTIONAL RESULTS
    diamond_files = glob.glob(os.path.join(DIAMOND_DIR, "*_hits_sial.tsv"))
    func_results = []

    for f in diamond_files:
        mid = os.path.basename(f).replace("_hits_sial.tsv", "")
        if os.path.getsize(f) > 0:
            hits = pd.read_csv(f, sep="\t", header=None)
            func_results.append({
                'mag_id': mid,
                'has_nanH': hits.iloc[:, 1].str.contains('nanH', case=False, na=False).any(),
                'has_nanE': hits.iloc[:, 1].str.contains('nanE', case=False, na=False).any(),
                'has_nanK': hits.iloc[:, 1].str.contains('nanK', case=False, na=False).any(),
                'has_nanA': hits.iloc[:, 1].str.contains('nanA', case=False, na=False).any(),
                'has_nanT': hits.iloc[:, 1].str.contains('nanT', case=False, na=False).any()
            })
        else:
            func_results.append({
                'mag_id': mid,
                'has_nanH': False, 'has_nanE': False, 'has_nanK': False, 'has_nanA': False, 'has_nanT': False
            })
            
    df_func = pd.DataFrame(func_results)
    
    for gene in GENES:
        df_func[gene] = df_func[gene].map({True: 'YES', False: 'NO'})
        
    print(f"Functional annotation loaded: {len(df_func)} MAGs")

    # 3. LOAD DESEQ2 RESULTS 
    if os.path.exists(DESEQ2_FILE):
        df_deseq = pd.read_csv(DESEQ2_FILE, sep="\t")
        
        def classify_enrichment(row):
            if pd.isna(row['padj']): return "NOT_SIGNIFICANT"
            if row['padj'] >= 0.05: return "NOT_SIGNIFICANT"
            if row['log2FoldChange'] > 0: return "ENRICHED_IN_AD"
            if row['log2FoldChange'] < 0: return "ENRICHED_IN_HC"
            return "NOT_SIGNIFICANT"
            
        df_deseq['deseq2_status'] = df_deseq.apply(classify_enrichment, axis=1)
        df_deseq = df_deseq[['Species', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'deseq2_status']]
        print(f"DESeq2 results loaded: {len(df_deseq)} Species")
    else:
        print("WARNING: DESeq2 file not found.")
        df_deseq = pd.DataFrame(columns=['Species', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'deseq2_status'])

    # 4. MERGE ALL TABLES
    df_final = pd.merge(df_tax, df_func, on='mag_id', how='left')
    
    # 4.2 Une com os resultados do DESeq2 usando a Species
    df_final = pd.merge(df_final, df_deseq, on='Species', how='inner') 

    # 5. CLEAN MISSING VALUES
    for gene in GENES:
        df_final[gene] = df_final[gene].fillna('NO')
    df_final['deseq2_status'] = df_final['deseq2_status'].fillna('NOT_TESTED')

    # 6. COLUMN ORDER 
    final_cols = [
        'mag_id', 'Phylum', 'Genus', 'Species',
        'has_nanH', 'has_nanE', 'has_nanK', 'has_nanA', 'has_nanT',
        'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'deseq2_status'
    ]
    final_cols = [c for c in final_cols if c in df_final.columns]
    df_final = df_final[final_cols]

    # 7. SAVE
    df_final.to_csv(OUTPUT_FILE, index=False, sep="\t")

    # 8. SUMMARY
    print(f"\nMaster table saved:\n{OUTPUT_FILE}")
    print(f"Total MAGs in final table: {len(df_final)}")
    print(f"Total Unique Species in final table: {df_final['Species'].nunique()}")
# RUN
if __name__ == "__main__":
    main()