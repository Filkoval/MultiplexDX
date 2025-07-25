import os
import pandas as pd
import numpy as np
import gseapy as gp
from gseapy.parser import read_gmt
from pathlib import Path
import re
import gseapy.plot
import hashlib

# Monkey-patch to truncate long filenames
def safe_filename(name, max_len=100):
    if len(name) > max_len:
        hash_part = hashlib.md5(name.encode()).hexdigest()
        return name[:max_len-9] + "_" + hash_part[:8]
    return name

orig_gseaplot = gseapy.plot.gseaplot

def patched_gseaplot(*args, **kwargs):
    ofname = kwargs.get("ofname")
    if ofname:
        base, ext = os.path.splitext(ofname)
        base = safe_filename(base)
        kwargs["ofname"] = f"{base}{ext}"
    return orig_gseaplot(*args, **kwargs)

gseapy.plot.gseaplot = patched_gseaplot


# === CONFIG ===
deg_path = r"C:\Users\sara.hrabovska\OneDrive - MultiplexDX, s.r.o\Dokumenty\Code\AIpredict\GSEA Visium\MDX_3007_histochoice\3007 histochoice kmeans4 DEGs GEX.xlsx"
sheet_name = "differential_expression sel col"

# Extract sample name from the DEG file folder (one level up from the file)
sample_name = os.path.basename(os.path.dirname(deg_path))


# Read Excel file, forcing first column (gene names) to be string type
df = pd.read_excel(deg_path, sheet_name=sheet_name, dtype={0: str})

# Load desired gene sets
gene_sets = {
    "Hallmark": "C:/Users/sara.hrabovska/OneDrive - MultiplexDX, s.r.o/Dokumenty/Code/AIpredict/GSEA Visium/MDX_3007_histochoice/msigdb_v2025.1.Hs_GMTs/h.all.v2025.1.Hs.symbols.gmt",
    "GO_BP":    "C:/Users/sara.hrabovska/OneDrive - MultiplexDX, s.r.o/Dokumenty/Code/AIpredict/GSEA Visium/MDX_3007_histochoice/msigdb_v2025.1.Hs_GMTs/c5.go.bp.v2025.1.Hs.symbols.gmt",
    "GO_MF":    "C:/Users/sara.hrabovska/OneDrive - MultiplexDX, s.r.o/Dokumenty/Code/AIpredict/GSEA Visium/MDX_3007_histochoice/msigdb_v2025.1.Hs_GMTs/c5.go.mf.v2025.1.Hs.symbols.gmt",
    "GO_CC":    "C:/Users/sara.hrabovska/OneDrive - MultiplexDX, s.r.o/Dokumenty/Code/AIpredict/GSEA Visium/MDX_3007_histochoice/msigdb_v2025.1.Hs_GMTs/c5.go.cc.v2025.1.Hs.symbols.gmt",
    "C8":       "C:/Users/sara.hrabovska/OneDrive - MultiplexDX, s.r.o/Dokumenty/Code/AIpredict/GSEA Visium/MDX_3007_histochoice/msigdb_v2025.1.Hs_GMTs/c8.all.v2025.1.Hs.symbols.gmt",
}

gene_sets_loaded = {
    name: read_gmt(str(path))
    for name, path in gene_sets.items()
}

# Use folder with space as you wanted
gsea_outdir_base = "./AIpredict/GSEA Visium/MDX_3007_histochoice/results"
os.makedirs(gsea_outdir_base, exist_ok=True)

# Load DEG data
deg_df = pd.read_excel(deg_path, sheet_name=sheet_name)
deg_df.rename(columns={'Feature Name': 'Gene'}, inplace=True)

# Clean gene names globally before analysis
deg_df['Gene'] = deg_df['Gene'].str.strip().str.upper()

def sanitize_term(term, max_length=100):
    # Replace illegal characters and truncate
    term = re.sub(r"[^\w\-]", "_", term)
    return term[:max_length]

# === RUN GSEA PRERANK FOR EACH CLUSTER AND GENE SET ===
clusters = [1, 2, 3, 4]
for cluster in clusters:
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

    # Sort by rank descending for prerank input
    rnk = df_cluster[['Gene', 'rank']].sort_values(by='rank', ascending=False)

    # Print the rank table to terminal
    print("Top 20 genes with their ranks for Cluster 1:")
    print(rnk.head(20))

    # Check for duplicates in the preranked input
    duplicates = rnk[rnk.duplicated(subset='Gene', keep=False)]
    if not duplicates.empty:
        print(f"⚠️ Found {duplicates['Gene'].nunique()} duplicate gene IDs in the ranking for cluster {cluster}:")
        print(duplicates.sort_values('Gene'))
    else:
        print(f"No duplicates found in preranked gene list for cluster {cluster}.")

    for gs_name, gs in gene_sets.items():
        outdir = os.path.join(gsea_outdir_base, f"gsea_cluster{cluster}_{gs_name}")
        prerank_dir = os.path.join(outdir, "prerank")  # <-- This is where GSEApy writes plots
        os.makedirs(prerank_dir, exist_ok=True)       # <-- Ensure it exists

        print(f"Running GSEA preranked for cluster {cluster} with gene set {gs_name}...")

        gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets_loaded[gs_name],
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

        
# === SUMMARIZE ALL GSEA RESULTS INTO ONE CSV ===
for gs_name in gene_sets.keys():
    all_results = []
    for cluster in clusters:
        folder = os.path.join(gsea_outdir_base, f"gsea_cluster{cluster}_{gs_name}")
        result_file = os.path.join(folder, "gseapy.gene_set.prerank.report.csv")

        if not os.path.exists(result_file):
            print(f"Missing: {result_file}")
            continue

        df = pd.read_csv(result_file)
        df['Cluster'] = cluster
        df['GeneSet'] = gs_name
        all_results.append(df)

    if all_results:
        summary_df = pd.concat(all_results, ignore_index=True)
        summary_csv_path = os.path.join(
            gsea_outdir_base,
            f"{sample_name}_GSEA_summary_{gs_name}_all_clusters.csv"
        )
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"✅ Summary saved for {gs_name} to: '{summary_csv_path}'")
    else:
        print(f"⚠️ No results found for gene set {gs_name}, skipping summary CSV.")