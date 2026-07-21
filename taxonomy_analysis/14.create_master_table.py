import pandas as pd
import os
import glob
import re

GTDB_DIR = "/home/marcos/PRJEB59406/table/gtdb_taxonomy/batch_*"
DIAMOND_DIR = "/home/marcos/PRJEB59406/table/mag_annotation/batch_*/diamond_matches"
DESEQ2_FILE = "/home/marcos/PRJEB59406/table/mag_abundance/deseq2_all_results.tsv"
OUTPUT_FILE = "/home/marcos/PRJEB59406/table/master_table_mag_sialidase_gtdb_deseq2.tsv"

# FUNCTIONS
def parse_taxonomy(tax_string):
    """
    Splits GTDB taxonomy into Phylum, Genus and Species
    """

    if pd.isna(tax_string):
        return pd.Series(['Unknown', 'Unknown', 'Unknown'])

    parts = tax_string.split(';')

    p = next((x.split('__')[1] for x in parts if 'p__' in x), 'Unknown')
    g = next((x.split('__')[1] for x in parts if 'g__' in x), 'Unknown')
    s = next((x.split('__')[1] for x in parts if 's__' in x), 'Unknown')

    return pd.Series([p, g, s])

# MAIN

def main():

    tax_files = glob.glob(
        os.path.join(GTDB_DIR, "gtdbtk.*.summary.tsv")
    )
    if not tax_files:
        print("ERROR: No GTDB-Tk files found.")
        return
    df_tax = pd.concat([
        pd.read_csv(
            f,
            sep="\t",
            usecols=['user_genome', 'classification']
        )
        for f in tax_files
    ])

    df_tax['mag_id'] = (
        df_tax['user_genome']
        .str.replace(".fa", "", regex=False)
    )
    # Parse taxonomy
    df_tax[['Phylum', 'Genus', 'Species']] = (
        df_tax['classification']
        .apply(parse_taxonomy)
    )
    print(f"Taxonomy loaded: {len(df_tax)} MAGs")

    # 2. LOAD DIAMOND FUNCTIONAL RESULTS

    diamond_files = glob.glob(
        os.path.join(DIAMOND_DIR, "*_hits_sial.tsv")
    )

    func_results = []

    for f in diamond_files:

        basename = os.path.basename(f)
        mid = basename.replace("_hits_sial.tsv", "")
        if os.path.getsize(f) > 0:
            hits = pd.read_csv(
                f,
                sep="\t",
                header=None
            )
            func_results.append({
                'mag_id': mid,
                'has_nanH':
                    'YES' if any(
                        hits.iloc[:, 1].str.contains(
                            'nanH',
                            case=False,
                            na=False
                        )
                    ) else 'NO',
                'has_nanE':
                    'YES' if any(
                        hits.iloc[:, 1].str.contains(
                            'nanE',
                            case=False,
                            na=False
                        )
                    ) else 'NO',
                'has_nanK':
                    'YES' if any(
                        hits.iloc[:, 1].str.contains(
                            'nanK',
                            case=False,
                            na=False
                        )
                    ) else 'NO',
                'has_nanA':
                    'YES' if any(
                        hits.iloc[:, 1].str.contains(
                            'nanA',
                            case=False,
                            na=False
                        )
                    ) else 'NO',
                'has_nanT':
                    'YES' if any(
                        hits.iloc[:, 1].str.contains(
                            'nanT',
                            case=False,
                            na=False
                        )
                    ) else 'NO',
                'gene_count': len(hits),
                'best_hit': hits.iloc[0, 1]
            })
        else:
            func_results.append({
                'mag_id': mid,
                'has_nanH': 'NO',
                'has_nanE': 'NO',
                'has_nanK': 'NO',
                'has_nanA': 'NO',
                'has_nanT': 'NO',
                'gene_count': 0,
                'best_hit': '-'
            })
    df_func = (
        pd.DataFrame(func_results)
        .drop_duplicates(subset='mag_id')
    )
    print(f"Functional annotation loaded: {len(df_func)} MAGs")

    # 3. LOAD DESEQ2 RESULTS

    if os.path.exists(DESEQ2_FILE):

        df_deseq = pd.read_csv(
            DESEQ2_FILE,
            sep="\t"
        )
        df_deseq['mag_id'] = df_deseq['MAG']

        def classify_enrichment(row):
            if pd.isna(row['padj']):
                return "NOT_SIGNIFICANT"
            if row['padj'] >= 0.05:
                return "NOT_SIGNIFICANT"
            if row['log2FoldChange'] > 0:
                return "ENRICHED_IN_AD"
            if row['log2FoldChange'] < 0:
                return "ENRICHED_IN_HC"
            return "NOT_SIGNIFICANT"
        df_deseq['deseq2_status'] = (
            df_deseq.apply(classify_enrichment, axis=1)
        )
        df_deseq = df_deseq[[
            'mag_id',
            'baseMean',
            'log2FoldChange',
            'lfcSE',
            'stat',
            'pvalue',
            'padj',
            'deseq2_status'
        ]]

        print(f"DESeq2 results loaded: {len(df_deseq)} MAGs")

    else:

        print("WARNING: DESeq2 file not found.")

        df_deseq = pd.DataFrame(columns=[
            'mag_id',
            'baseMean',
            'log2FoldChange',
            'lfcSE',
            'stat',
            'pvalue',
            'padj',
            'deseq2_status'
        ])
    # 4. MERGE ALL TABLES

    df_final = pd.merge(
        df_tax,
        df_func,
        on='mag_id',
        how='left'
    )
    df_final = pd.merge(
        df_final,
        df_deseq,
        on='mag_id',
        how='left'
    )
    # 5. CLEAN MISSING VALUES

    for gene in [
        'has_nanH',
        'has_nanE',
        'has_nanK',
        'has_nanA',
        'has_nanT'
    ]:
        df_final[gene] = df_final[gene].fillna('NO')
    df_final['gene_count'] = (
        df_final['gene_count']
        .fillna(0)
    )
    df_final['best_hit'] = (
        df_final['best_hit']
        .fillna('-')
    )
    df_final['deseq2_status'] = (
        df_final['deseq2_status']
        .fillna('NOT_TESTED')
    )

    # 6. COLUMN ORDER

    final_cols = [
        'mag_id',
        'Phylum',
        'Genus',
        'Species',
        'has_nanH',
        'has_nanE',
        'has_nanK',
        'has_nanA',
        'has_nanT',
        'gene_count',
        'best_hit',
        'baseMean',
        'log2FoldChange',
        'lfcSE',
        'stat',
        'pvalue',
        'padj',
        'deseq2_status'
    ]
    df_final = df_final[final_cols]

    # 7. SAVE

    df_final.to_csv(
        OUTPUT_FILE,
        index=False,
        sep="\t"
    )
    # 8. SUMMARY


    print(f"Master table saved:")
    print(OUTPUT_FILE)
    print("\nGene summary:")
    print(f"nanH-positive MAGs: {len(df_final[df_final['has_nanH'] == 'YES'])}")
    print(f"nanE-positive MAGs: {len(df_final[df_final['has_nanE'] == 'YES'])}")
    print(f"nanK-positive MAGs: {len(df_final[df_final['has_nanK'] == 'YES'])}")
    print(f"nanA-positive MAGs: {len(df_final[df_final['has_nanA'] == 'YES'])}")
    print(f"nanT-positive MAGs: {len(df_final[df_final['has_nanT'] == 'YES'])}")
    print("\nDESeq2 summary:")
    print(f"Enriched in AD: {len(df_final[df_final['deseq2_status'] == 'ENRICHED_IN_AD'])}")
    print(f"Enriched in HC: {len(df_final[df_final['deseq2_status'] == 'ENRICHED_IN_HC'])}")
# RUN

if __name__ == "__main__":
    main()