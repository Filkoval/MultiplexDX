# ====== LOAD PACKAGES AND ARGUMENTS =======
cat("Starting script...\n")
flush.console()

args <- commandArgs(trailingOnly = TRUE)
required_pkgs <- c("patchwork", "Seurat", "ggplot2","RColorBrewer","jsonlite","magick")
invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

DIR = "C:\\MultiplexDX-local\\Presentation_ecotyper\\Space_ranger-3007x"
FILE = "C:\\MultiplexDX-local\\Presentation_ecotyper\\Space_ranger-3007x\\analysis\\clustering\\gene_expression_kmeans_4_clusters\\clusters.csv"
CLUSTER = 3
OUTPUT_TIS = 'C:\\MultiplexDX-local\\Presentation_ecotyper\\plots\\3007x\\sample_3007x_cluster_1\\3007x_tissue_1.pdf'
OUTPUT_NO = 'C:\\MultiplexDX-local\\Presentation_ecotyper\\plots\\3007x\\sample_3007x_cluster_1\\3007x_no_tissue_1.pdf'
OUTPUT_ANNOT = 'C:\\MultiplexDX-local\\Presentation_ecotyper\\plots\\3007x\\sample_3007x_cluster_1\\3007x_annot_1.png'
POINT = 6

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

#new_img <- Read10X_Image(
#  image.dir = "path/to/new/image",
#  filter.matrix = TRUE
#)

# Priradiť ho namiesto starého
#seurat_obj@images[["slice1"]] <- new_img

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

ggsave(OUTPUT_TIS, plot = plot1)
ggsave(OUTPUT_NO, plot = plot2)




