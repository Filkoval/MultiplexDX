import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# === CONFIG ===
deg_path = r"C:\Users\sara.hrabovska\OneDrive - MultiplexDX, s.r.o\Dokumenty\Code\AIpredict\GSEA Visium\MDX_3007_histochoice\3007 histochoice kmeans4 DEGs GEX.xlsx"
clusters = [1, 2, 3, 4]
gene_set_names = ["Hallmark", "GO_BP", "GO_MF", "GO_CC", "C8"]

# Extract folder and sample name
deg_folder = os.path.dirname(deg_path)
sample_name = os.path.basename(deg_folder)
gsea_outdir_base = os.path.join(deg_folder, "results")  # <- updated path to match subfolder location

# === GENERATE SEPARATE PDFs FOR EACH GENE SET ===
for gs_name in gene_set_names:
    summary_csv_path = os.path.join(gsea_outdir_base, f"{sample_name}_GSEA_summary_{gs_name}_all_clusters.csv")

    if not os.path.exists(summary_csv_path):
        print(f"⚠️ Summary CSV not found for {gs_name}, skipping: {summary_csv_path}")
        continue
    else:
        print(f"✅ Found: {summary_csv_path}")

    # Load gene set specific summary CSV
    df = pd.read_csv(summary_csv_path)

    # Validate columns exist
    required_cols = {'Term', 'NES', 'NOM p-val', 'Cluster'}
    if not required_cols.issubset(df.columns):
        print(f"⚠️ Missing columns in {gs_name} summary CSV, skipping.")
        continue

    # Prepare output PDF path for this gene set
    pdf_path = os.path.join(gsea_outdir_base, f"{sample_name}_GSEA_{gs_name}.pdf")

    with PdfPages(pdf_path) as pdf:

        # Pivot to create NES and p-value matrices
        nes_df = df.pivot(index='Term', columns='Cluster', values='NES')
        pval_df = df.pivot(index='Term', columns='Cluster', values='NOM p-val')

        # Filter: keep terms significant in at least one cluster
        sig_terms = pval_df.lt(0.05).any(axis=1)
        nes_df = nes_df.loc[sig_terms]
        pval_df = pval_df.loc[sig_terms]

        # Replace infs with NaNs
        nes_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        pval_df.replace([np.inf, -np.inf], np.nan, inplace=True)

        # Drop rows with NaNs
        nes_df.dropna(how='any', inplace=True)
        pval_df = pval_df.loc[nes_df.index]

        if nes_df.empty:
            print(f"⚠️ No valid NES data for {gs_name}, skipping.")
            continue

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
            continue

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
        plt.suptitle(f"{gs_name} Gene Set Enrichment (NES, p < 0.05)", y=1.02, fontsize=10)
        pdf.savefig(g.fig, bbox_inches='tight')
        plt.close(g.fig)

    print(f"✅ PDF saved to: {pdf_path}")
