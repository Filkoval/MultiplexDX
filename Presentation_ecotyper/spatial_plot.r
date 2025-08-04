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

#CROP PICTURE
crop <- function(
)
{
  img <- image_read(paste0(DIR,"\\annotation.png"))
  path =  paste0(DIR,"/spatial/tissue_positions.csv")
  csv_file <- read.csv(path, header = TRUE, sep = ",")
  
  up <- min(csv_file$array_row, na.rm = TRUE)
  left <- min(csv_file$array_col, na.rm = TRUE)
  down <- max(csv_file$array_row, na.rm = TRUE)
  right <- max(csv_file$array_col, na.rm = TRUE)

  print(down)

  up_left <- (csv_file[csv_file$array_row == up & csv_file$array_col == left, ])
  down_right <- (csv_file[csv_file$array_row == down & csv_file$array_col == right, ])


  info <- image_info(image_read(paste0(DIR,"\\annotation.png")))
  scale_x <- info$width / 11350
  scale_y <- info$height / 19290

  x=21097*scale_x
  y=20564*scale_y

  # Definuj rohy (napr. ľavý horný = 50,50; pravý dolný = 250,350)
  x1 <- (down_right$pxl_col_in_fullres)*scale_x
  y1 <-  (down_right$pxl_row_in_fullres)*scale_y
  x2 <- (up_left$pxl_col_in_fullres)*scale_x
  y2 <- (up_left$pxl_row_in_fullres)*scale_y

  #3007 xylene5_Scan1 (4,w=11350, h=19290
  
  

  # Vypočítaj šírku, výšku a posun
  w <- round(x2 - x1 + 1)
  h <- round(y2 - y1 + 1)
  ox <- round(x1-x)
  oy <- round(y1-y)


  # Orez v magick (geometria: "wxh+ox+oy")
  img <- image_crop(img,geometry = sprintf("%dx%d+%d+%d", w, h, ox, oy))
  print(img)
  image_write(img,, path = paste0(OUTPUT_ANNOT))
}
crop()


