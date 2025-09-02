#SEURAT SETUP
 'This short script shows the basic functionalities of the Seurat library for spatial transcriptomics analyses.'

#INPUT PARAMETERS
args <- commandArgs(trailingOnly = TRUE)
  # args[1] - Path to output folder(Where the output of this vignette will be saved into new Seurat_output directory)
  # args[2] - Name of sample(Any, but no whitespaces!)
  # args[3] - Space Ranger output directory path

  # RUN BY: Rscript path_to_this_script args[1] "args[2]" args[3]
  # e.g.(working when the working directory is the exact directory from example below) Rscript ./Seurat_vignette.r . "skuska" ./Cell_ranger_outputs/sample1

  # Example directory:
  # |--- Seurat_vignette.r
  # |--- Seurat_output
  # |     |--- plots(here save the output plots, does not have to be crated beforehand)
  # |     |--- rds_objects(here save temporary files which fasten data processing in the future runs, does not have to be crated beforehand)
  # |     |--- cropped_output.png(cropped image by fiducials, does not have to be crated beforehand)
  # |--- Cell_Ranger_outputs
  # |     |--- sample1(is \outs directory created by space ranger from sample1)


if(length(args) < 3) { #Testing regime
  print("Testing regime on.")
  PATH <- ".\\Seurat_output"
  NAME <- "Sample1"
  INPUT <- ".\\Cell_ranger_outputs\\sample1"
} else {
  print("Arguments read.")
  PATH <- paste0(args[1],"\\Seurat_output")
  NAME <- args[2]
  INPUT <- args[3]
}
if(!dir.exists(PATH)){dir.create(PATH)}
if(!dir.exists(paste0(PATH,"/rds_objects"))){dir.create(paste0(PATH,"/rds_objects"))}
if(!dir.exists(paste0(PATH,"/plots"))){dir.create(paste0(PATH,"/plots"))}

#LIBRARIES

print("Loading libraries...")
required_pkgs <- c("jsonlite","patchwork", "Seurat", "ggplot2", "dplyr", "glmGamPoi", "rlang", "grid","imagefx","magick") #EBImage needed to be installed by BiocManager
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE) 
}

# Načítanie balíkov so skrytím štartovacích správ
invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
print("All libraries loaded")

#DATA LOADING
load_data <- function(
)
{
  print("Loading data...")
  sample <- Load10X_Spatial(
    data.dir = INPUT, #creates Seurat object
  )
  sample[["orig.ident"]] <- "sample1"
  sample@project.name <- "sample1"

  sample <- SetAssayData( #Adds default count for SCT
    sample,
    assay    = "Spatial",                        
    layer    = "counts",                     
    new.data = sample@assays$Spatial@layers[["counts.Gene Expression"]]
  )

  print("Data loaded")
  return(sample)
}

#RAW UMI COUNT VISUALIZATION
show_data <- function(
  object #Seurat object
)
{
    print("Creating plots...")
    plot1 <- VlnPlot(object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
    plot2 <- SpatialFeaturePlot(object, features = "nCount_Spatial") + theme(legend.position = "right")
    plots <- wrap_plots(plot1, plot2)
    p(plots,"raw_counts")
}

#COMPARISON OF SCTRANFORMATION AND NORMALIZATION
comparison_of_normalization <- function(
  object #Seurat object
)
{
    print("Launching comparison...")
    temp <- transform(object,"brain") #SCT

    print("Normalizing...") 
    temp <- NormalizeData(object, verbose = FALSE, assay = "Spatial") #Classic normalization

    print("Correlating groups...")
    temp <- GroupCorrelation(temp, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
    temp <- GroupCorrelation(temp, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

    print("Creating plots...")
    p1 <- GroupCorrelationPlot(brain, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
    theme(plot.title = element_text(hjust = 0.5))
    p2 <- GroupCorrelationPlot(brain, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
    theme(plot.title = element_text(hjust = 0.5))
    p(wrap_plots(p1, p2),"comparison_SCT_and_norm")
}

#CREATES GRADIENT PLOT - VISUALIZATION OF DESIRED GENES
plot_grad <- function(
  transformed, #Seurat object with SCT Assay
  num = 5, #number of most frequented genes to plot
  inter = FALSE, #interactivity(true if interactive environemnt needed)
  genes = FALSE #FALSE for most frequented genes plotting mode, character vector for specific genes plotting mode
)
{
  print("Creating plot...")
   
  if(!genes){genes <- most_expressed(transformed, num)}
  if(is.character(genes))
  {
    plot <- SpatialFeaturePlot(transformed, features = genes, interactive = inter, pt.size.factor = 4.2)
    return(plot)
  }
  else{print("Invalid input in genes, should be character vector - c().")}
}

#CREATES GROUP PLOT - VISUALIZATION OF GENES GROUPED BY METADATA COLUMNS
plot_group <- function(
  transformed, #Seurat object with SCT Assay
  group = "orig.ident", #metadata column, to group spots by(see available by 'print(colnames(sample1@meta.data))')
  inter = FALSE #TRUE for interactive mode
  )
{
  print("Creating plot...")
   
  plot <- SpatialDimPlot(transformed, group.by = group, interactive = inter)
  p(plot,"clustering_by_group_by")
}

#SCTRANSFORMATION - IMPORTANT FIRST STEP FOR OTHER PROCESES
transform <- function(
  object, #Seurat object to transform
  name = "Temporal_rds_file" #name to save transformation preprocessed data by
)
{
  file_name = paste0(PATH,"/rds_objects/",name,"_transformed")
  if(!file.exists(file_name))
  {
    print("SCTransforming...")
    temp <- SCTransform(object, assay = "Spatial", verbose = FALSE)
    saveRDS(temp, file = file_name)
    print("SCT assay created, rds file created.")
    return(temp)
  }
  else
  {
    print("Transformed data already available.")
    temp <- readRDS(file_name)
    return(temp)
  }
}

#PCA - DATA TRANSFORMATION TO FIND GENE DEPENDENCIES
create_PCA <-function(
  object #Seurat object
)
{
  file_name = paste0(PATH,"/rds_objects/",NAME,"_PCA")
  if(!file.exists(file_name))
  {
    print("Running PCA...")
    temp <- RunPCA(object, assay = "SCT", verbose = FALSE)
    saveRDS(temp, file = file_name)
    print("PCA done, rds file created.")
  }
  else
  {
    print("PCA already available.")
    temp <- readRDS(file_name)
  }
  return(temp)
}

#GRAPH CLUSTERING AND VISUALIZING
create_clusters <- function(
  object #Seurat object
)
{
  file_name = paste0(PATH,"/rds_objects/",NAME,"_Clusters")
  if(!file.exists(file_name))
  {
    print("Clustering...")
    temp <- create_PCA(object) #run PCA
    temp <- FindNeighbors(temp, reduction = "pca", dims = 1:30) #connects dependent genes and creates a graph with valued edges, based on shared neigbours
    #(dims selects number of PCA - less = not enoug dependency, too much = misleads to false dependencies)
    temp <- FindClusters(temp, verbose = FALSE) #finds clusters - Louvain-Leiden greedy algoritmus
    temp <- RunUMAP(temp, reduction = "pca", dims = 1:30) #2D visualization of gene dependencies
    saveRDS(temp, file = file_name)
  }
  else
  {
    print("Clusters already available.")
    temp <- readRDS(file_name)
  }
  return(temp)
}

#VISUALIZE CLUSTERS
visualize_clusters <- function(
  object #Seurat object
  )
{
   
  temp <- object
  p1 <- DimPlot(temp, reduction = "umap", label = TRUE) #visualizes umap reduction - clustering of gene dependencies
  p2 <- SpatialDimPlot(temp, label = TRUE, label.size = 3) #automatically visualizes all seurat_clusters vector done by create_clusters
  p3 <- SpatialDimPlot(temp, cells.highlight = CellsByIdentities(object = temp, idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3) #viusalizes desired clusters(2,1,4,3,5,8)
  p(p1,name="visualized_clusters_1")
  p(p2,name="visualized_clusters_2") 
  p(p3,name="visualized_clusters_3") 
}

#INTERACTIVE SELECTION
interactive_selection_graph <- function(
  object #Seurat object
  )
{
  temp <- create_clusters(object) #creates clusters
  p <- LinkedDimPlot(temp) #interactive graph with spot selection
  print(p)
}

#DIFFERENTIAL EXPRESSION
diff_exp <- function(
  object, #Seurat object with clustering - stored in Idents
  cluster1, #number of first cluster
  cluster2 #number of second cluster
)
{
   
  temp <- object
  de_markers <- FindMarkers(temp, ident.1 = cluster1, ident.2 = cluster2) #finds and sorts differentially expressed genes between the two clusters
  plot <- SpatialFeaturePlot(
  object = temp,
  features = rownames(de_markers)[1:3],
  alpha = c(0.1, 1),
  ncol = 3
) & theme(
  legend.position = "bottom",
  legend.key.width  = unit(0.2, "cm"),
  legend.key.height = unit(0.3, "cm"),
  legend.text       = element_text(size = 6),
  legend.title      = element_text(size = 7)
)
}

#VISUALIZATION OF GENERAL PLOTS
p <- function(
  plot,
  name = "temporal_picture"
)
{
  dev.new()
  print(plot)
  ggsave(paste0(PATH, "/plots/", name, ".pdf"), plot = plot)
}

#COEXPRESSION
coexpress <- function(
  object, #Seurat object 
  gene_list #list of genes to coexpress
  ) {
  temp <- AddModuleScore(object, features = gene_list, name = "myModule")
  plot1 <- SpatialFeaturePlot(temp, features = "myModule1", pt.size.factor = 2.7)
  plot2 <- SpatialFeaturePlot(object = object, features = unlist(gene_list), pt.size.factor = 2.7)

  my_genes <- list()
  for (gene in gene_list) {
    my_genes[[gene]] <- object[["SCT"]]@counts[gene, ]
  }

  expr_mat <- do.call(rbind, my_genes)

  nonzero_spots <- colSums(expr_mat != 0) == nrow(expr_mat)
  spots_num <- sum(nonzero_spots)

  plot <- plot1 / plot2 +
    plot_annotation(
      title = paste0("Počet spotov, kde je expresia všetkých génov významná: ", spots_num)
    )
  p(plot, "coexpression")
}


#CROP PICTURE
crop <- function(
)
{
  img <- image_read(paste0(INPUT,"\\spatial\\tissue_hires_image.png"))
  path =  paste0(INPUT,"/spatial/tissue_positions.csv")
  csv_file <- read.csv(path, header = TRUE, sep = ",")
  scale_factors <- fromJSON(paste0(INPUT,"\\spatial\\scalefactors_json.json"))
  scale <- scale_factors$tissue_hires_scalef
  
  

  up <- min(csv_file$array_row, na.rm = TRUE)
  left <- min(csv_file$array_col, na.rm = TRUE)
  down <- max(csv_file$array_row, na.rm = TRUE)
  right <- max(csv_file$array_col, na.rm = TRUE)

   print(down)

  up_left <- (csv_file[csv_file$array_row == up & csv_file$array_col == left, ])
  down_right <- (csv_file[csv_file$array_row == down & csv_file$array_col == right, ])

  # Definuj rohy (napr. ľavý horný = 50,50; pravý dolný = 250,350)
  x1 <- (down_right$pxl_col_in_fullres)*scale
  y1 <-  (down_right$pxl_row_in_fullres)*scale
  x2 <- (up_left$pxl_col_in_fullres)*scale
  y2 <- (up_left$pxl_row_in_fullres)*scale

  # Vypočítaj šírku, výšku a posun
  w <- round(x2 - x1 + 1)
  h <- round(y2 - y1 + 1)
  ox <- round(x1)
  oy <- round(y1)

  # Orez v magick (geometria: "wxh+ox+oy")
  img <- image_crop(img,geometry = sprintf("%dx%d+%d+%d", w, h, ox, oy))
  print(img)
  image_write(img,, path = paste0(PATH,"/cropped_output.png"))
}

#MOST EXPRESSED GENES
most_expressed <- function(
  object, #Seurat object
  num #creates list of "num" most expressed genes
)
{
  avg_expr <- rowMeans(sample1[["SCT"]]@counts)
  genes <- names(sort(avg_expr, decreasing = TRUE)[1:num])
  print(genes)
  return(genes)
}

#MULTIPLE SLICES
multiple <- function(
  slice1, #First slice to merge
  slice2 #SEcond slice to merge
)
{
  merged <- merge(slice1, slice2)
  DefaultAssay(merged) <- "SCT"
  VariableFeatures(merged) <- c(VariableFeatures(slice1), VariableFeatures(slice2))
  merged <- RunPCA(merged, verbose = FALSE)
  merged <- FindNeighbors(merged, dims = 1:30)
  merged <- FindClusters(merged, verbose = FALSE)
  merged <- RunUMAP(merged, dims = 1:30)
  p(DimPlot(merged, reduction = "umap", group.by = c("ident", "orig.ident")))
}


#COMMENT & UNCOMMENT FUNCTIONS BY ADDING/REMOVING #

sample1 <- load_data()
sample1 <- transform(sample1,NAME)
#show_data(sample1)
#p(plot_grad(sample1,inter=TRUE))
#sample1 <- create_clusters(sample1)
#visualize_clusters(sample1)
#diff_exp(sample1,2,4)
#interactive_selection_graph(Sample1)
coexpress(sample1,sample1[["SCT"]]@var.features[1:2])
crop()
multiple(sample1, sample1)
print("End of script.")