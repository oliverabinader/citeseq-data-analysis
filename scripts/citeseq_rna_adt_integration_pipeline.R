############################################################
## CITE-seq (RNA + ADT) Full Processing Pipeline
############################################################

library(Seurat)
library(cowplot)
library(dplyr)
library(dsb)
library(ggplot2)

############################################################
## Load all samples (each folder = one sample)
############################################################

seurat_objects <- list()

for (folder in list.files()) {

  sample_data <- Read10X(data.dir = file.path(folder, "outs/filtered_feature_bc_matrix"))

  seurat_obj <- CreateSeuratObject(
    counts = sample_data$`Gene Expression`,
    project = folder,
    min.cells = 3
  )

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj[["ADT"]] <- CreateAssay5Object(counts = sample_data$`Antibody Capture`)

  seurat_objects[[folder]] <- seurat_obj
}

############################################################
## QC filtering (kept as project-specific thresholds)
############################################################

seurat_objects <- lapply(seurat_objects, function(obj) {
  subset(
    obj,
    subset =
      nFeature_RNA >= 200 &
      nFeature_RNA <= 6000 &
      percent.mt <= 5 &
      nCount_RNA <= 50000 &
      nCount_ADT <= 10000
  )
})

############################################################
## Split into NM (nasal mucosa) and Lung datasets
############################################################
## NOTE: folder names define biological grouping

seurat_objects_NM <- seurat_objects[c("Folder_C1", "Folder_C3", "Folder_C9", "Folder_C11")]
seurat_objects_Lung <- seurat_objects[c("Folder_C2", "Folder_C4", "Folder_C10", "Folder_C12")]

############################################################
## RNA integration (SCT) - NM
############################################################

options(future.globals.maxSize = 10000 * 1024^2)

for (i in 1:length(seurat_objects_NM)) {
  seurat_objects_NM[[i]] <- SCTransform(seurat_objects_NM[[i]], verbose = FALSE)
}

nm.features <- SelectIntegrationFeatures(seurat_objects_NM, nfeatures = 3000)
seurat_objects_NM <- PrepSCTIntegration(seurat_objects_NM, anchor.features = nm.features, verbose = FALSE)

nm.anchors <- FindIntegrationAnchors(
  object.list = seurat_objects_NM,
  normalization.method = "SCT",
  anchor.features = nm.features,
  verbose = FALSE
)

nm.integrated <- IntegrateData(
  anchorset = nm.anchors,
  normalization.method = "SCT",
  verbose = FALSE,
  new.assay.name = "SCTintegrated"
)

nm.integrated <- RunPCA(nm.integrated, reduction.name = "pca_SCTintegrated", reduction.key = "pca_", assay = "SCTintegrated", verbose = FALSE)
#ElbowPlot(nm.integrated, ndims = 30, reduction = "pca_SCTintegrated")
nm.integrated <- RunUMAP(nm.integrated, dims = 1:30, reduction = "pca_SCTintegrated", reduction.name = "umap_SCTintegrated", reduction.key = "umap_", assay = "SCTintegrated")

nm.integrated.plot <- DimPlot(nm.integrated, group.by = "orig.ident", reduction = "umap_SCTintegrated") & theme(legend.position = "right") 

############################################################
## RNA integration (SCT) - Lung
############################################################

for (i in 1:length(seurat_objects_Lung)) {
  seurat_objects_Lung[[i]] <- SCTransform(seurat_objects_Lung[[i]], verbose = FALSE)
}

lung.features <- SelectIntegrationFeatures(seurat_objects_Lung, nfeatures = 3000)
seurat_objects_Lung <- PrepSCTIntegration(seurat_objects_Lung, anchor.features = lung.features, verbose = FALSE)

lung.anchors <- FindIntegrationAnchors(
  object.list = seurat_objects_Lung,
  normalization.method = "SCT",
  anchor.features = lung.features
)

lung.integrated <- IntegrateData(
  anchorset = lung.anchors,
  normalization.method = "SCT",
  new.assay.name = "SCTintegrated",
  verbose = FALSE
)

lung.integrated <- RunPCA(lung.integrated, reduction.name = "pca_SCTintegrated", reduction.key = "pca_", assay = "SCTintegrated", verbose = FALSE)
lung.integrated <- RunUMAP(lung.integrated, dims = 1:30, reduction = "pca_SCTintegrated", reduction.name = "umap_SCTintegrated", reduction.key = "umap_", assay = "SCTintegrated")

lung.integrated.plot <- DimPlot(lung.integrated, group.by = "orig.ident", reduction = "umap_SCTintegrated") & theme(legend.position = "right") 

############################################################
## ADT processing - NM
############################################################

seurat_objects_NM <- lapply(seurat_objects_NM, function(obj) {

  DefaultAssay(obj) <- "ADT"

  VariableFeatures(obj) <- rownames(obj[["ADT"]]) # Use all ADT features for dimensional reduction

  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE) %>%
    ScaleData(features = VariableFeatures(obj)) %>%
    RunPCA(features = VariableFeatures(obj))

  obj
})

#ElbowPlot(object = seurat_objects_NM$Folder_C1, ndims = 19)

nm.adt.anchors <- FindIntegrationAnchors(
  object.list = seurat_objects_NM,
  reduction = "rpca",
  dims = 1:19
)

nm.adt.integrated <- IntegrateData(
  anchorset = nm.adt.anchors,
  dims = 1:19,
  new.assay.name = "ADTintegrated"
)

nm.adt.integrated <- ScaleData(nm.adt.integrated, verbose = FALSE)
#for very large number of protein antibodies (for a large panel), add the following parameters: ScaleData(do.scale = FALSE,  do.center = TRUE)

nm.adt.integrated <- RunPCA(nm.adt.integrated, reduction.name = "pca_ADTintegrated", verbose = FALSE, reduction.key = "pca_", assay = "ADTintegrated")
nm.adt.integrated <- RunUMAP(nm.adt.integrated, dims = 1:19, reduction = "pca_ADTintegrated", reduction.name = "umap_ADTintegrated", reduction.key = "umap_", assay = "ADTintegrated")

nm.adt.integrated.plot <- DimPlot(nm.adt.integrated, group.by = "orig.ident", reduction = "umap_ADTintegrated") & theme(legend.position = "right") 

############################################################
## ADT processing - Lung
############################################################

seurat_objects_Lung <- lapply(seurat_objects_Lung, function(obj) {

  DefaultAssay(obj) <- "ADT"

  VariableFeatures(obj) <- rownames(obj[["ADT"]])

  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE) %>%
    ScaleData(features = VariableFeatures(obj)) %>%
    RunPCA(features = VariableFeatures(obj))

  obj
})

lung.adt.anchors <- FindIntegrationAnchors(
  object.list = seurat_objects_Lung,
  reduction = "rpca",
  dims = 1:19
)

lung.adt.integrated <- IntegrateData(
  anchorset = lung.adt.anchors,
  dims = 1:19,
  new.assay.name = "ADTintegrated"
)

lung.adt.integrated <- ScaleData(lung.adt.integrated, verbose = FALSE)

lung.adt.integrated <- RunPCA(lung.adt.integrated, reduction.name = "pca_ADTintegrated", verbose = FALSE, reduction.key = "pca_", assay = "ADTintegrated")
lung.adt.integrated <- RunUMAP(lung.adt.integrated, dims = 1:19, reduction = "pca_ADTintegrated", reduction.name = "umap_ADTintegrated", reduction.key = "umap_", assay = "ADTintegrated")

lung.adt.integrated.plot <- DimPlot(lung.adt.integrated, group.by = "orig.ident", reduction = "umap_ADTintegrated") & theme(legend.position = "right") 

############################################################
## Add ADT integration into RNA-integrated objects
############################################################

nm.integrated[["ADTintegrated"]] <- nm.adt.integrated[["ADTintegrated"]]
nm.integrated[["pca_ADTintegrated"]] <- nm.adt.integrated[["pca_ADTintegrated"]]
nm.integrated[["umap_ADTintegrated"]] <- nm.adt.integrated[["umap_ADTintegrated"]]

lung.integrated[["ADTintegrated"]] <- lung.adt.integrated[["ADTintegrated"]]
lung.integrated[["pca_ADTintegrated"]] <- lung.adt.integrated[["pca_ADTintegrated"]]
lung.integrated[["umap_ADTintegrated"]] <- lung.adt.integrated[["umap_ADTintegrated"]]

#################################################################
## Saving Seurat object before performing WNN (for future use)
#################################################################

saveRDS(nm.integrated, "/path/to/WNN_Analysis/nm.integrated.pre_wnn.rds")
saveRDS(lung.integrated, "/path/to/WNN_Analysis/lung.integrated.pre_wnn.rds")

############################################################
## Multimodal neighbor graph (Weighted_Nearest_Neighbor)
############################################################

find_multimodal_neighbors <- function(seurat_obj) {
  seurat_obj <- FindMultiModalNeighbors(
    seurat_obj, 
    reduction.list = list("pca_SCTintegrated", "pca_ADTintegrated"), 
    dims.list = list(1:30, 1:19), 
    modality.weight.name = c("SCTintegrated.weight", "ADTintegrated.weight")
  )
  return(seurat_obj)
}

nm.wnn <- find_multimodal_neighbors(nm.integrated)
lung.wnn <- find_multimodal_neighbors(lung.integrated)

############################################################
## Downstream clustering + UMAP
############################################################

#Define a function to run UMAP, find clusters, and generate a plot

run_UMAP_find_clusters <- function(obj) {
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = FALSE)
  plot <- DimPlot(obj, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 3.5)
  return(list(obj = obj, plot = plot))
}

nm.results <- run_UMAP_find_clusters(nm.wnn)
lung.results <- run_UMAP_find_clusters(lung.wnn)

nm.wnn <- nm.results$obj
nm.wnn.plot <- nm.results$plot

lung.wnn <- lung.results$obj
lung.wnn.plot <- lung.results$plot

############################################################
## SingleR annotation (Lung samples example shown)
############################################################

library(celldex)
library(SingleR)

ref <- ImmGenData() #this reference has some Lung samples compatible/relevant for this study

lung_counts <- GetAssayData(lung.wnn, layer = "counts", assay = "SCT") #matrix containing counts expression data (raw data) for use as the "test" dataset in SingleR
 
pred <- SingleR(test = lung_counts, ref = ref, assay.type.test = 1, labels = ref$label.main) #assay.type.test = 1 means use the count data not normalized one

#Adding cell types to our Seurat object
lung.wnn@meta.data$singleR_label <- pred$labels

#Cell Type visualization
ggsave(
  filename = "wnn_umap.Lung.label.tiff", 
  plot = DimPlot(lung.wnn, reduction = "wnn.umap", group.by = "singleR_label", label = TRUE, repel = TRUE, label.size = 3.5) + ggtitle(""), 
  device = "tiff", 
  path = "/path/to/celltype_identification/", 
  width = 8,
  height = 6,
  units = "in",
  dpi = 350,
  compression = "lzw"
)

#################################################################
## Saving Seurat object before performing WNN (for future use)
#################################################################

saveRDS(nm.wnn, "/path/to/WNN_Analysis/nm.integrated.after_wnn.rds")
saveRDS(lung.wnn, "/path/to/WNN_Analysis/lung.integrated.after_wnn.rds")

############################################################
## Differential Expression Analysis (For all clusters)
############################################################

# Prepare object to run differential analysis 
nm.wnn <- PrepSCTFindMarkers(nm.wnn)
lung.wnn <- PrepSCTFindMarkers(lung.wnn)

# Define clustering identity
Idents(nm.wnn) <- "wsnn_res.0.3"
Idents(lung.wnn) <- "wsnn_res.0.3"

#Run differential expression
nm_markers <- FindAllMarkers(nm.wnn)
lung_markers <- FindAllMarkers(lung.wnn)

############################################################
## Generate cloupe file from Seurat object
############################################################

library(loupeR)

#Please note that only the count matrix (RNA assay) will be used to generate the CLOUPE file. SC transformations, integrated or other assays will not be stored in the Loupe file.
#Reference: https://www.10xgenomics.com/support/software/loupe-browser/latest/tutorials/introduction/lb-louper

create_loupe_from_seurat(
  nm.wnn, 
  output_name = "loupe_from_NM_samples", 
  output_dir = "/path/to/cloupe/", 
  force = TRUE
)

create_loupe_from_seurat(
  lung.wnn, 
  output_name = "loupe_from_Lung_samples", 
  output_dir = "/path/to/cloupe/", 
  force = TRUE
)

############################################################
## ADT DotPlot Visualization (NM)
############################################################

DefaultAssay(nm.wnn) <- "ADTintegrated"

nm_sample_titles <- c(
  "Folder_1" = "WT_0_NM",
  "Folder_3" = "KO_0_NM",
  "Folder_9" = "KO_14_NM",
  "Folder_11" = "WT_14_NM"
)

for (sample in names(nm_sample_titles)) {
  p <- DotPlot(
    subset(nm.wnn, subset = orig.ident == sample),
    features = rownames(nm.wnn),
    cols = c("blue", "red"),
    dot.scale = 8
  ) +
    RotatedAxis() +
    ggtitle(nm_sample_titles[sample])

  ggsave(
    filename = paste0("DotPlot_", sample, "_NM.tiff"),
    plot = p,
    path = "/path/to/dotplot/ADT/",
    device = "tiff",
    width = 8,
    height = 6,
    units = "in",
    dpi = 350,
    compression = "lzw"
  )
}

############################################################
## ADT DotPlot (Lung)
############################################################

DefaultAssay(lung.wnn) <- "ADTintegrated"

lung_sample_titles <- c(
  "Folder_2" = "WT_0_Lung",
  "Folder_4" = "KO_0_Lung",
  "Folder_10" = "KO_14_Lung",
  "Folder_12" = "WT_14_Lung"
)

for (sample in names(lung_sample_titles)) {
  p <- DotPlot(
    subset(lung.wnn, subset = orig.ident == sample),
    features = rownames(lung.wnn),
    cols = c("blue", "red"),
    dot.scale = 8
  ) +
    RotatedAxis() +
    ggtitle(lung_sample_titles[sample])

  ggsave(
    filename = paste0("DotPlot_", sample, "_Lung.tiff"),
    plot = p,
    path = "/path/to/dotplot/ADT/",
    device = "tiff",
    width = 8,
    height = 6,
    units = "in",
    dpi = 350,
    compression = "lzw"
  )
}
