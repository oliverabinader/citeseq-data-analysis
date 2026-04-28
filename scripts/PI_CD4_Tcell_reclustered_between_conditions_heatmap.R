#####################################################################################################################
## PI-driven Analysis: T cell Cluster Heatmap (WT vs KO)
#####################################################################################################################
## Description:
## This script performs a custom downstream analysis based on PI-provided cell barcode annotations for T cell subsets
## (e.g., CD4 or CD8 T cells).
##
## The input CSV file contains:
## - Cell barcodes
## - Predefined cluster assignments
## - Experimental condition labels (e.g., WT/B6 vs KO/BLT)
##
## Workflow:
## 1. Subset Seurat object using provided barcodes
## 2. Map external cluster/condition annotations to cells
## 3. Extract gene expression from SCT-integrated assay
## 4. Compute average gene expression per cluster (pseudo-bulk aggregation)
## 5. Perform gene-wise Z-score scaling
## 6. Generate a heatmap using ComplexHeatmap with hierarchical clustering (row and column-wise) across clusters
##
## Output:
## - Heatmap comparing cluster-level expression profiles between WT and KO conditions
##
## Notes:
## - Clusters are NOT derived from Seurat clustering but are externally defined by the PI
## - This analysis enables direct comparison of matched clusters (e.g., Cluster 1 WT vs Cluster 1 KO)
## - Can be adapted for CD4 or CD8 T cell subsets by changing the input annotation file
#####################################################################################################################

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

############################################################
## Load Seurat object (post-WNN recommended)
############################################################

lung_obj <- readRDS("/path/to/lung.integrated.after_wnn.rds") 

############################################################
## Load PI-provided barcode annotation file
############################################################

#  Example: CD4 T cell clusters
barcode_metadata <- read.csv(
  "/path/to/CD4_Tcell_clusters_PI_annotations.csv"
)

############################################################
## Subset Seurat object to relevant barcodes
############################################################

# It’s good to verify which barcodes are present in our data

valid_barcodes <- barcode_metadata$Barcode[
  barcode_metadata$Barcode %in% rownames(lung_obj@meta.data)
]
# Sanity check: length(valid_barcodes) == length(barcode_metadata$Barcode)

lung_subset <- subset(lung_obj, cells = valid_barcodes) 

############################################################
## Add PI cluster annotations to metadata
############################################################

# Before moving forward, you may need to split the reclustered column in the annotation file into say 
# condition and cluster, so that you could add these seperatly into the metadata.

# Then, you assign condition and cluster to cell barcode in seurat object. Here, we are just adding clusters.
# cluster_column containing cluster + condition labels
lung_subset@meta.data$PI_cluster <- barcode_metadata[[cluster_column]][
  match(rownames(lung_subset@meta.data), barcode_metadata$Barcode)
]

############################################################
## Extract expression matrix (SCT integrated assay)
############################################################

DefaultAssay(lung_subset) <- "SCTintegrated"

exp_matrix <- as.matrix(lung_subset[["SCTintegrated"]]$data) #if you encounter here an error, you may need to read the data in the form of a dataframe.
meta <- lung_subset@meta.data

############################################################
## Compute average expression per cluster
############################################################

# Get the unique cluster labels (as they appear in your matrix colnames)
# cluster_labels <- unique(gsub("\\.\\d+$", "", colnames(exp_matrix)))  # remove trailing ".1", ".2", etc if any.
# Note: colnames of exp_matrix may contain same colnames as meta$Condition_ClusteringfromPI 
cluster_labels <- unique(colnames(exp_matrix))

cluster_avg_list <- list()

for (label in cluster_labels) {
  matching_cols <- which(colnames(exp_matrix) == label)
  cluster_avg <- rowMeans(exp_matrix[, matching_cols, drop = FALSE])
  cluster_avg_list[[label]] <- cluster_avg
}

avg_df <- as.data.frame(cluster_avg_list)
#avg_df$Gene <- rownames(exp_matrix)
#avg_df <- avg_df[, c("Gene", setdiff(names(avg_df), "Gene"))]  # Move Gene to first column
#rownames(avg_df) <- avg_df$Gene
#avg_df$Gene <- NULL  # Optional

############################################################
## Scale data (Z-score per gene)
############################################################

scaled_data <- t(scale(t(as.matrix(avg_df)), center = TRUE, scale = TRUE))

############################################################
## Generate Heatmap
############################################################

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#col_fun(seq(-3, 3))

scaled_data <- t(scale(t(as.matrix(avg_df)), center = TRUE, scale = TRUE))
heatmap_plot <- Heatmap(
  matrix = scaled_data,
  col = col_fun,
  column_title = "PI-defined clusters",
  row_title = "Genes",
  name = "Expression",
  show_row_dend = FALSE,
  show_row_names = FALSE,
  use_raster = TRUE,
  #row_dend_width = unit(10, "cm"),
  #clustering_distance_rows = "pearson", 
  #clustering_method_rows = "single"
)

############################################################
## Save Heatmap
############################################################

pdf("/path/to/output/PI_Tcell_cluster_heatmap.pdf", width = 10, height = 8)
draw(heatmap_plot)
dev.off()
