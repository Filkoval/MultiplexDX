# ====== LOAD PACKAGES AND ARGUMENTS =======
cat("Starting script...\n")
flush.console()

args <- commandArgs(trailingOnly = TRUE)
required_pkgs <- c("patchwork", "Seurat", "ggplot2","RColorBrewer","jsonlite","magick")
invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))



DIR = args[1] 
FILE = args[2]
CLUSTER = args[3]
OUTPUT_TIS = args[4]
OUTPUT_NO = args[5]
POINT = args[6]




# ======= LOAD OBJECT ========
seurat_obj <- Load10X_Spatial(data.dir = DIR)
clust <- read.csv(FILE)

rownames(clust) <- clust$Barcode
seurat_obj$kmeans <- clust[Cells(seurat_obj), "cluster"]
Idents(seurat_obj) <- "kmeans"
seurat_obj@meta.data$kmeans <- clust[colnames(seurat_obj), "Cluster"]
seurat_obj@meta.data$highlight <- ifelse(seurat_obj@meta.data$kmeans == CLUSTER,
                                         paste0("Cluster-", CLUSTER),
                                         "Other")
seurat_obj@meta.data$highlight <- factor(seurat_obj@meta.data$highlight, levels = c(paste0("Cluster-", CLUSTER), "Other"))

Idents(seurat_obj) <- "highlight"
cluster_name<-paste0("Cluster-", CLUSTER)

# ========= PLOT ========
cols <- c(
   cluster_name= "#440024ff",
   "Other"= "#8DA0CB"                   
)
names(cols) <- c(cluster_name,"Other")

plot1 <- SpatialDimPlot(
  seurat_obj,
  label = FALSE, 
  pt.size.factor = as.numeric(POINT)-1, 
  cols = cols, 
  image.alpha = 0.85, 
  cols.highlight = c("#DE2D26", "grey50")
  ) + NoLegend()


plot2 <- SpatialDimPlot(
seurat_obj,
label = FALSE, 
pt.size.factor = as.numeric(POINT), 
cols = cols, 
image.alpha = 0, 
cols.highlight = c("#DE2D26", "grey50"),
) + NoLegend()

# =========== UMAP =============
seurat_obj <- SetAssayData( #Adds default count for SCT
    seurat_obj,
    assay    = "Spatial",                        
    layer    = "counts",                     
    new.data = seurat_obj@assays$Spatial@layers[["counts.Gene Expression"]]
  )

umap_df <- read.csv(paste0(DIR,"/analysis/umap/gene_expression_2_components/projection.csv"))
kmeans_df <- read.csv(FILE)

merged_df <- merge(umap_df, kmeans_df, by.x = "Barcode", by.y = "Barcode")
merged_df <- merged_df[match(colnames(seurat_obj), merged_df$Barcode), ]
seurat_obj$kmeans_cluster <- as.factor(merged_df$Cluster)
umap_mat <- as.matrix(merged_df[, c("UMAP.1", "UMAP.2")])

rownames(umap_mat) <- merged_df$Barcode
umap_mat <- umap_mat[colnames(seurat_obj), ]

seurat_obj[["umap_spaceranger"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key = "UMAPSR_",
  assay = DefaultAssay(seurat_obj)
)

seurat_obj$kmeans_cluster <- as.factor(merged_df$Cluster)
plotik <- DimPlot(seurat_obj, reduction = "umap_spaceranger", group.by = "kmeans_cluster", label = TRUE) +
  ggtitle("SpaceRanger UMAP s K-means klastrami")

ggsave(OUTPUT_TIS, plot = plot1)
ggsave(OUTPUT_NO, plot = plot2)
ggsave(paste0(DIR, "/UMAP.pdf"), plot = plotik)
print("konec")

