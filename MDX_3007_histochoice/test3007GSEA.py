import os
import sys
import pandas as pd
import numpy as np
import gseapy as gp
from gseapy.parser import read_gmt
from pathlib import Path
import re
import gseapy.plot
import hashlib
import yaml

# === CONFIG ===================================================
# first argument - ptah to YAML config file in following format:
# -----------------------------------------
#"xlsx_url": <PATH>
#"sheet_name": <NAME>
#"gene_set_names" : <PYTHON LIST OF GENES>
#"urls": <PYTHON LIST OF PATHS>
#"output_dir" : <PATH>
#"clusters" : <NUMBER OF CLUSTERS>
# -----------------------------------------
# Read data
INPUT_FILE = sys.argv[1]
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(BASE_DIR,"../"))
print(str(os.getcwd()))
with open(INPUT_FILE, 'r', encoding='utf-8') as f:
    data = yaml.safe_load(f)
# ==============================================================

deg_path = data["xlsx_url"]
sheet_name = data["sheet_name"]

sample_name = os.path.basename(os.path.dirname(deg_path)) # Extract sample name from the DEG file folder (one level up from the file)
df = pd.read_excel(deg_path, sheet_name=sheet_name, dtype={0: str}) # Read Excel file, forcing first column (gene names) to be string type

def correct_filename(filename, max_len=200):
    base, ext = os.path.splitext(filename)
    if len(base) > max_len:
        hash_part = hashlib.md5(base.encode()).hexdigest()
        base = base[:max_len-9] + "_" + hash_part[:8]
    filename = f"{base}{ext}"
    return filename

def load_data():
    # Load desired gene sets
    gene_sets = {}
    for gene,url in zip(data["gene_set_names"],data["urls"]):
        gene_sets[gene] = read_gmt(str(url))


    # Use folder with space as you wanted
    gsea_outdir_base = data["output_dir"]
    os.makedirs(gsea_outdir_base, exist_ok=True)

    # Load DEG data
    deg_df = pd.read_excel(deg_path)
    deg_df.rename(columns={'Feature Name': 'Gene'}, inplace=True)
    return gene_sets, gsea_outdir_base, deg_df

def create_cluster_dataframe(cluster):
    log2fc_col = f'Cluster {cluster} Log2 fold change'
    padj_col = f'Cluster {cluster} Adjusted p value'
    
    # Use full dataframe (no filtering)
    df_cluster = deg_df[['Gene', log2fc_col, padj_col]].dropna()

    # Replace zeros in adjusted p-values with a small number (e.g. 1e-300)
    df_cluster[padj_col] = df_cluster[padj_col].clip(lower=1e-300)

       # Calculate ranking
    df_cluster['rank'] = -np.log10(df_cluster[padj_col]) * np.sign(df_cluster[log2fc_col])

    # Add tiny random jitter to break ties in rank values
    np.random.seed(42)  # for reproducibility
    df_cluster['rank'] += np.random.uniform(-1e-6, 1e-6, size=len(df_cluster))
    
    # Remove duplicates keeping highest absolute rank
    df_cluster['abs_rank'] = df_cluster['rank'].abs()
    df_cluster = df_cluster.sort_values(by='abs_rank', ascending=False)
    df_cluster = df_cluster.drop_duplicates(subset='Gene', keep='first')
    df_cluster = df_cluster.drop(columns='abs_rank')
    return df_cluster

def control_duplicates(duplicates):
    if not duplicates.empty:
        print(f"⚠️ Found {duplicates['Gene'].nunique()} duplicate gene IDs in the ranking for cluster {cluster}:")
        print(duplicates.sort_values('Gene'))
    else:
        print(f"No duplicates found in preranked gene list for cluster {cluster}.")

def top_genes():
    # Print the rank table to terminal
    print("Top 20 genes with their ranks for Cluster 1:")
    print(rnk.head(20))

def create_output_dir(cluster,gs_name):
    outdir = os.path.join(gsea_outdir_base, f"gsea_cluster{cluster}_{gs_name}")
    prerank_dir = os.path.join(outdir, "prerank")  # <-- This is where GSEApy writes plots
    os.makedirs(prerank_dir, exist_ok=True)
    return outdir

def prerank():

    print(f"Running GSEA preranked for cluster {cluster} with gene set {gs_name}...")

    gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets[gs_name],
            outdir=Path(outdir).as_posix(),
            no_plot=True,
            min_size=15,
            max_size=500,
            permutation_num=100,
            seed=42,
            verbose=True,
            format='png',  # PNG avoids long filenames and PDF issues
            format_short=True
        )
    
def result_file(cluster,gs_name):
    folder = os.path.join(gsea_outdir_base, f"gsea_cluster{cluster}_{gs_name}")
    result_file = os.path.join(folder, "gseapy.gene_set.prerank.report.csv")

    if not os.path.exists(result_file):
        print(f"Missing: {result_file}")
        return None

    df = pd.read_csv(result_file)
    df['Cluster'] = cluster
    df['GeneSet'] = gs_name
    
    return df

def create_summary():
    summary_df = pd.concat(all_results, ignore_index=True)
    summary_csv_path = os.path.join(
        gsea_outdir_base,
        f"{sample_name}_GSEA_summary_{gs_name}_all_clusters.csv"
    )
    summary_df.to_csv(summary_csv_path, index=False)
    print(f"✅ Summary saved for {gs_name} to: '{summary_csv_path}'")
    

gene_sets, gsea_outdir_base, deg_df = load_data()
# Clean gene names globally before analysis
deg_df['Gene'] = deg_df['Gene'].str.strip().str.upper()


# === RUN GSEA PRERANK FOR EACH CLUSTER AND GENE SET ===
for i in range(data["clusters"]):

    cluster = i+1
    df_cluster = create_cluster_dataframe(cluster)
    # Sort by rank descending for prerank input
    rnk = df_cluster[['Gene', 'rank']].sort_values(by='rank', ascending=False)

    top_genes()
    # Check for duplicates in the preranked input
    duplicates = rnk[rnk.duplicated(subset='Gene', keep=False)]
    control_duplicates(duplicates)

    for gs_name, gs in gene_sets.items():

        outdir = create_output_dir(cluster,gs_name)
        prerank()

        
# === SUMMARIZE ALL GSEA RESULTS INTO ONE CSV ===
for gs_name in gene_sets.keys():
    all_results = []
    for i in range(data["clusters"]):
        cluster = i+1
        all_results.append(result_file(cluster,gs_name))

    if all_results:
        create_summary()
    else:
        print(f"⚠️ No results found for gene set {gs_name}, skipping summary CSV.")