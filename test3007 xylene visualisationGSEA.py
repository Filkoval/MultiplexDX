import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# === CONFIG ===
deg_path = r"C:\Users\sara.hrabovska\OneDrive - MultiplexDX, s.r.o\Dokumenty\Code\AIpredict\GSEA Visium\MDX_3007_xylene\3007 xylene kmeans4 DEGs GEX.xlsx"
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
        fdr_df = df.pivot(index='Term', columns='Cluster', values='FDR q-val')
        sig_terms = fdr_df.lt(0.25).any(axis=1)  # Use 0.25 or 0.05 depending on stringency
        nes_df = nes_df.loc[sig_terms]
        pval_df = pval_df.loc[sig_terms]
        fdr_df = fdr_df.loc[sig_terms]

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
        annot_df = fdr_df.map(lambda x: "*" if x < 0.25 else "")

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
        plt.suptitle(f"{gs_name} Gene Set Enrichment (NES, filtered by FDR q-value < 0.25)", y=1.02, fontsize=10)
        pdf.savefig(g.fig, bbox_inches='tight')
        plt.close(g.fig)

    print(f"✅ PDF saved to: {pdf_path}")


    # Convert data to long-form for plotting
long_df = df[df['Term'].isin(nes_df.index)].copy()

# Only keep significant ones by FDR
long_df = long_df[long_df['FDR q-val'] < 0.25]

# Calculate -log10(FDR) for size
long_df['logFDR'] = -np.log10(long_df['FDR q-val'].clip(lower=1e-10))  # avoid -inf
long_df['NES_clipped'] = long_df['NES'].clip(-3, 3)  # for color scaling

# Set up the figure
plt.figure(figsize=(12, max(6, 0.4 * long_df['Term'].nunique())))
sns.set(style="whitegrid")

# Create the bubble plot
bubble_plot = sns.scatterplot(
    data=long_df,
    x="Cluster",
    y="Term",
    size="logFDR",
    hue="NES_clipped",
    palette="coolwarm",
    sizes=(50, 500),
    edgecolor="black",
    linewidth=0.3,
    legend='brief'
)

# Adjust aesthetics
bubble_plot.set_title(f"{gs_name} Gene Set Enrichment (NES color, bubble size ~ FDR)", fontsize=14)
bubble_plot.tick_params(axis='y', labelsize=8)
bubble_plot.tick_params(axis='x', labelsize=10)
bubble_plot.set_xlabel("Cluster", fontsize=12)
bubble_plot.set_ylabel("Gene Set Term", fontsize=12)

# Improve font clarity
for label in bubble_plot.get_yticklabels():
    label.set_fontsize(10)
    label.set_fontweight('bold')

# Move legend outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title="NES / -log10(FDR)")

# Tight layout and save
plt.tight_layout()
pdf.savefig(bubble_plot.get_figure(), bbox_inches='tight')
plt.close()

