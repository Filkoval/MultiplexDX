import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import yaml


# === CONFIG ===================================================
# first argument - ptah to YAML config file in following format:
# -----------------------------------------
# "xlsx path": <PATH>
# "clusters": <NUMBER OF CLUSTERS>
# "gene_set_names": <PYTHON LIST OF NAMES>
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
clusters = [(i+1) for i in range(int(data["clusters"]))]
gene_set_names = data["gene_set_names"]

# Extract folder and sample name
deg_folder = os.path.dirname(deg_path)
sample_name = os.path.basename(deg_folder)
gsea_outdir_base = os.path.join(deg_folder, "results")  # <- updated path to match subfolder location


def load_gene_set(): #load gene set from file
    fpath = f"{sample_name}_GSEA_summary_{gs_name}_all_clusters.csv"
    summary_csv_path = os.path.join(gsea_outdir_base, fpath)
    print(summary_csv_path)
    if not os.path.exists(summary_csv_path):
        print(f"⚠️ Summary CSV not found for {gs_name}, skipping: {summary_csv_path}")
        return None
    else:
        print(f"✅ Found: {summary_csv_path}")
        df = pd.read_csv(summary_csv_path)
        return df


def column_exist(): # Validate columns exist
    required_cols = {'Term', 'NES', 'NOM p-val', 'Cluster'}
    if not required_cols.issubset(df.columns):
        print(f"⚠️ Missing columns in {gs_name} summary CSV, skipping.")
        return False
    return True


def Reshape():
        
    # Pivot to create NES and p-value matrices
    nes_df = df.pivot(index='Term', columns='Cluster', values='NES') # enrichment score
    pval_df = df.pivot(index='Term', columns='Cluster', values='NOM p-val') # possibility of random statistical closenss

    # Filter: keep terms significant in at least one cluster
    sig_terms = pval_df.lt(0.05).any(axis=1) #filteras out stochaistic genes
    nes_df = nes_df.loc[sig_terms] # filtering rows based on True/False series above 
    pval_df = pval_df.loc[sig_terms]

    # Replace infs with NaNs
    nes_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    pval_df.replace([np.inf, -np.inf], np.nan, inplace=True)

    # Drop rows with NaNs
    nes_df.dropna(how='any', inplace=True)
    pval_df = pval_df.loc[nes_df.index]

    if nes_df.empty:
        print(f"⚠️ No valid NES data for {gs_name}, skipping.")
        return None
    return nes_df,pval_df


def define_clusters(nes_df,pval_df):

    # Ensure columns (clusters) are in the correct order
    for c in clusters:
        if c not in nes_df.columns:
            nes_df[c] = np.nan
            pval_df[c] = np.nan

    nes_df = nes_df[clusters]
    pval_df = pval_df[clusters]

    # Recreate annotation stars AFTER cleanup
    annot_df = pval_df.map(lambda x: "*" if x < 0.05 else "")

    # Check for enough data for clustering
    if nes_df.shape[0] < 2 or nes_df.shape[1] < 2:
        print(f"⚠️ Not enough data to cluster for {gs_name} (shape: {nes_df.shape}), skipping.")
        return None
    return annot_df


def graph_design():
    sns.set(font_scale=0.5)
    g = sns.clustermap(
        nes_df,
        cmap='coolwarm',
        center=0,
        annot=annot_df,
        fmt="",
        linewidths=0.3,
        col_cluster=True,
        row_cluster=True,
        method='average',
        metric='correlation',
        cbar_kws={'label': 'NES'},
        figsize=(10, max(6, 0.25 * len(nes_df)))  # Adjust size based on term count
    )
    return g


def draw_graphs():

    plt.suptitle(f"{gs_name} Gene Set Enrichment (NES, p < 0.05)", y=1.02, fontsize=10)
    pdf.savefig(graph_design().figure, bbox_inches='tight')
    plt.close(graph_design().figure)

# === GENERATE SEPARATE PDFs FOR EACH GENE SET ===
for gs_name in gene_set_names:

    # Prepare output PDF path for this gene set
    pdf_path = os.path.join(gsea_outdir_base, f"{sample_name}_GSEA_{gs_name}.pdf")
    df = load_gene_set()
    if df is None: continue
    if not column_exist(): continue

    with PdfPages(pdf_path) as pdf:
        val = Reshape()
        if val is None : continue
        nes_df,pval_df = val
        annot_df = define_clusters(nes_df,pval_df)
        if annot_df is None : continue
        draw_graphs()

    print(f"✅ PDF saved to: {pdf_path}")
