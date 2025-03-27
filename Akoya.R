rm(list=ls())

# Due to its large size, the input file "Measurements_P21.txt" could not be uploaded to this repository. 
# The file is deposited in the KU Leuven Research Data Repository (https://doi.org/10.48804/992X8C) and it is available upon request.

# Setting wd and projectfolder for saving plots
set.seed(2)
projectFolder <- "img_01/"
dir.create(projectFolder, showWarnings = FALSE)

# Loading packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(tidyverse)
library(RColorBrewer)


###### Part 1: Process segmented Akoya data ###### 

# Load Akoya segmentation data
HNSCC_Akoya_Measures <- read_tsv("Measurements_P21.txt") 

# Display dimensions of the loaded data
dim(HNSCC_Akoya_Measures) #1,764,150    1,242

# Create empty data frame for metadata
HNSCC_Akoya_meta <- data.frame(matrix(ncol = 0, nrow = nrow(HNSCC_Akoya_Measures)))
HNSCC_Akoya_meta$CellIDs <- paste0("CellID_",rownames(HNSCC_Akoya_meta))
HNSCC_Akoya_meta$Image <- "P21"

# Extract mean intensities of each marker for each cell
HNSCC_Akoya_Measures_filter <- HNSCC_Akoya_Measures[,grepl("Cell: Mean", colnames(HNSCC_Akoya_Measures))]
colnames(HNSCC_Akoya_Measures_filter) <- gsub(": Cell: Mean","",colnames(HNSCC_Akoya_Measures_filter))

# Extract metadata (first 8 columns)
HNSCC_Akoya_meta2 <- HNSCC_Akoya_Measures[,c(1:8)]
HNSCC_Akoya_meta <- cbind(HNSCC_Akoya_meta2, HNSCC_Akoya_meta)

# Display dimensions of processed data
dim(HNSCC_Akoya_Measures_filter) #1,764,150      61
dim(HNSCC_Akoya_meta) #1,764,150      10

# Convert marker data to a data frame
HNSCC_Akoya_Measures_filter <- as.data.frame(HNSCC_Akoya_Measures_filter)

# Assign row names for consistency
rownames(HNSCC_Akoya_Measures_filter) <- HNSCC_Akoya_meta$CellIDs
rownames(HNSCC_Akoya_meta) <- HNSCC_Akoya_meta$CellIDs

# Not using DAPI as a gene expression value, yet as a meta.data value
# DAPI is a measure of cell/nucleus-size, so it might affect final clustering if implemented in the expression matrix 
HNSCC_Akoya_meta$DAPI <- HNSCC_Akoya_Measures_filter$DAPI
HNSCC_Akoya_Measures_filter$DAPI <- NULL

# Create Seurat object
Akoya_Seurat_test <- CreateSeuratObject(counts = t(HNSCC_Akoya_Measures_filter), 
                                        project = "Akoya_Scan1", 
                                        assay = "Akoya_Mean", 
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = HNSCC_Akoya_meta)

Akoya_Seurat_test
# An object of class Seurat
# 60 features across 1764150 samples within 1 assay
# Active assay: Akoya_Mean (60 features, 0 variable features)
#  1 layer present: counts

# Plot XY scatter of all cells
png(paste0(projectFolder, "/codex_XY_scatter_P21.png"), height=1200, width=2000)
print(FeatureScatter(Akoya_Seurat_test, 
               feature1 = "Centroid X µm", 
               feature2 = "Centroid Y µm", 
               raster = FALSE, 
               shuffle = TRUE, 
               group.by = "Name",
               pt.size = 0.1))
dev.off()


# Note: DAPI has been removed as an expression feature for downstream steps, DAPI intensities can be found in the meta.data
Akoya_Features <- rownames(Akoya_Seurat_test[["Akoya_Mean"]])
Akoya_Features
#  [1] "CD44"            "CD20"            "CD107a"          "CD45RO"          "CD31"            "Collagen IV"     "Podoplanin"      "SMA"             "CD4"             "E-cadherin"      "Vimentin"        "CD11c"           "GranzymeB"       "EpCAM"
# [15] "CD3e"            "CD34"            "CD79a"           "CD68"            "CD141"           "FOXP3"           "MPO"             "Pan-Cytokeratin" "IDO1"            "CD8"             "PD-1"            "CD163"           "PD-L1"           "CD39"
# [29] "CD45"            "HLA-E"           "CD138"           "IFNG"            "HLA-DR"          "TIGIT"           "LAG3"            "VISTA"           "ICOS"            "b-catenin1"      "Beta-actin"      "Bcl-2"           "CD21"            "Caveolin"
# [43] "CD66"            "GATA3"           "Keratin8/18"     "HIF1A"           "SOX2"            "Ki67"            "HistoH3 Phospho" "T-bet"           "LIF"             "TFAM"            "TCF-1"           "TOX"             "PCNA"            "Granzyme K"
# [57] "Runx3"           "TP63"            "CX3CR1"          "CD56"

# Display number of markers
length(Akoya_Features) #60

# Check the intensities of all markers
pdf(paste0(projectFolder, "/codex_QC_intensities_P21.pdf"), width = 20, height = 30)
print(VlnPlot(Akoya_Seurat_test, features = Akoya_Features, pt.size = 0))
dev.off()

# Use default normalization
Akoya_Seurat_sample1 <- NormalizeData(Akoya_Seurat_test, normalization.method = "CLR", margin = 2)

# Regress for nCount_Akoya_mean (overal signal strength)
Akoya_Seurat_sample1 <- ScaleData(Akoya_Seurat_sample1, assay = "Akoya_Mean", vars.to.regress = c("nCount_Akoya_Mean"))




# Extract centroid coordinates from metadata
centroids <- Akoya_Seurat_sample1@meta.data[,c("Centroid Y µm", "Centroid X µm", "CellIDs")]

# Create centroid object
cents <- SeuratObject::CreateCentroids(centroids)

# Define segmentation data
segmentations_data <- list(
    "centroids" = cents,
    "segmentation" = NULL
  )

# Try creating FOV object with improved error handling
result_coords <- tryCatch(
  CreateFOV(
    coords = segmentations_data,
    type = c("segmentation", "centroids"),
    molecules = NULL,
    assay = "Spatial"
  ),
  warning = function(w) {
    message("Warning encountered: ", conditionMessage(w))
    NULL
  },
  error = function(e) {
    message("Error encountered: ", conditionMessage(e))
    NULL
  }
)

# Store result in Seurat object
Akoya_Seurat_sample1[["Akoya"]] <- result_coords





###### Part 2: Determine marker scores ###### 

# Define marker lists for TLS-like and Type-1 immunity cells
TLS_like_TME_markers <- c("Caveolin","CD11c","CD141","CD163","CD20","CD21","CD31","CD34","CD3e","CD4","CD40","CD68","CD79a","CX3CR1","IDO1","IFNG","Ki67","PCNA","PD-1","PD-L1")
Type1_immunity_markers <- c("Caveolin","CD11c","CD141","CD163","CD31","CD34","CD3e","CD4","CD66","CD68","CD8","FOXP3","GranzymeB","HLA-A","IFNG","Ki67","LAG3","PCNA","PD-1","PD-L1")

# Identify shared and unique markers
Shared_markers <- intersect(TLS_like_TME_markers, Type1_immunity_markers)
TLS_like_TME_markers_unique <- setdiff(TLS_like_TME_markers, Shared_markers)
Type1_immunity_markers_unique <- setdiff(Type1_immunity_markers, Shared_markers)

# Display unique markers
TLS_like_TME_markers_unique # "CD20"   "CD21"   "CD40"   "CD79a"  "CX3CR1" "IDO1"
Type1_immunity_markers_unique #"CD66"      "CD8"       "FOXP3"     "GranzymeB" "HLA-A"     "LAG3"

# Identify missing markers in this sample
setdiff(TLS_like_TME_markers, Akoya_Features) #CD40
setdiff(Type1_immunity_markers, Akoya_Features) #HLA-A

# Define hub scores based on marker intensities
Akoya_Seurat_sample1 <- AddhubScore(Akoya_Seurat_sample1, features = list(TLS_like_TME_markers), ctrl = 5, name = "TLS_like_TME_score", nbin = 10, slot = "scale.data")
Akoya_Seurat_sample1 <- AddhubScore(Akoya_Seurat_sample1, features = list(Type1_immunity_markers), ctrl = 5, name = "Type1_immunity_score", nbin = 10, slot = "scale.data")
Akoya_Seurat_sample1 <- AddhubScore(Akoya_Seurat_sample1, features = list(TLS_like_TME_markers_unique), ctrl = 5, name = "TLS_like_TME_markers_unique", nbin = 10, slot = "scale.data")
Akoya_Seurat_sample1 <- AddhubScore(Akoya_Seurat_sample1, features = list(Type1_immunity_markers_unique), ctrl = 5, name = "Type1_immunity_markers_unique", nbin = 10, slot = "scale.data")

# Summarize hub scores
summary(Akoya_Seurat_sample1@meta.data$TLS_like_TME_score1)
summary(Akoya_Seurat_sample1@meta.data$Type1_immunity_score1)
summary(Akoya_Seurat_sample1@meta.data$TLS_like_TME_markers_unique1)
summary(Akoya_Seurat_sample1@meta.data$Type1_immunity_markers_unique1)

# Optional: reverse y-coordinates for visualization
Akoya_Seurat_sample1_viz <- Akoya_Seurat_sample1
Akoya_Seurat_sample1_viz@meta.data$`Centroid Y µm` <- -Akoya_Seurat_sample1_viz@meta.data$`Centroid Y µm` #to change meta.data 
Akoya_Seurat_sample1_viz@images$Akoya@boundaries$centroids@coords[,1] <- -Akoya_Seurat_sample1_viz@images$Akoya@boundaries$centroids@coords[,1] #to change FOV 

# Visualizing each score separately
dir.create(paste0(projectFolder, "imagefeatureplots"), showWarnings = FALSE)

png(paste0(projectFolder, "imagefeatureplots/codex_sample_P21_TLS_like_TME_markers.png"), height=2000, width=3000)
print(ImageFeaturePlot(Akoya_Seurat_sample1_viz,"TLS_like_TME_score1", min.cutoff = "q10", max.cutoff = "q90", size = 1.5, alpha = 0.9, border.color = "NA", dark.background = F))
dev.off()

png(paste0(projectFolder, "imagefeatureplots/codex_sample_P21_Type1_immunity_markers.png"), height=2000, width=3000)
print(ImageFeaturePlot(Akoya_Seurat_sample1_viz,"Type1_immunity_score1", min.cutoff = "q10", max.cutoff = "q90", size = 1.5, alpha = 0.9, border.color = "NA", dark.background = F))
dev.off()

png(paste0(projectFolder, "imagefeatureplots/codex_sample_P21_TLS_like_TME_markers_unique.png"), height=2000, width=3000)
print(ImageFeaturePlot(Akoya_Seurat_sample1_viz,"TLS_like_TME_markers_unique1", min.cutoff = "q10", max.cutoff = "q90", size = 1.5, alpha = 0.9, border.color = "NA", dark.background = F))
dev.off()

png(paste0(projectFolder, "imagefeatureplots/codex_sample_P21_Type1_immunity_markers_unique.png"), height=2000, width=3000)
print(ImageFeaturePlot(Akoya_Seurat_sample1_viz,"Type1_immunity_markers_unique1", min.cutoff = "q10", max.cutoff = "q90", size = 1.5, alpha = 0.9, border.color = "NA", dark.background = F))
dev.off()

# Visualizing two scores blended
dir.create(paste0(projectFolder, "imagefeatureplots_blended"), showWarnings = FALSE)

png(paste0(projectFolder, "imagefeatureplots_blended/codex_sample_P21_blended.png"), height=2000, width=3000)
print(ImageFeaturePlot(Akoya_Seurat_sample1_viz, features = c("TLS_like_TME_score1","Type1_immunity_score1"), min.cutoff = "q10", max.cutoff = "q90", size = 1.5, alpha = 0.9, border.color = "NA", dark.background = F, blend = TRUE))
dev.off()

png(paste0(projectFolder, "imagefeatureplots_blended/codex_sample_P21_unique_blended.png"), height=2000, width=3000)
print(ImageFeaturePlot(Akoya_Seurat_sample1_viz, features = c("TLS_like_TME_markers_unique1","Type1_immunity_markers_unique1"), min.cutoff = "q10", max.cutoff = "q90", size = 1.5, alpha = 0.9, border.color = "NA", dark.background = F, blend = TRUE))
dev.off()




###### Part 3.1: Focus on TLS-like example section ###### 

# First ensure correct column names for centroids
colnames(Akoya_Seurat_sample1@meta.data)[10] <- "Centroid.X.µm"
colnames(Akoya_Seurat_sample1@meta.data)[11] <- "Centroid.Y.µm"

# Define the region of interest for TLS-like
Akoya_TLS_like_subset <- subset(Akoya_Seurat_sample1, Centroid.X.µm > 3500 & Centroid.X.µm < 4160 & Centroid.Y.µm > 11900 & Centroid.Y.µm < 12450)

Akoya_TLS_like_subset
# An object of class Seurat
# 60 features across 5087 samples within 1 assay
# Active assay: Akoya_Mean (60 features, 60 variable features)
#  3 layers present: counts, data, scale.data
#  1 spatial field of view present: Akoya

# Optional: reverse y-coordinates for visualization
Akoya_TLS_like_subset@meta.data$`Centroid.Y.µm` <- -Akoya_TLS_like_subset@meta.data$`Centroid.Y.µm` #to change meta.data 
Akoya_TLS_like_subset@images$Akoya@boundaries$centroids@coords[,1] <- -Akoya_TLS_like_subset@images$Akoya@boundaries$centroids@coords[,1] #to change FOV 

# Plot XY scatter of all cells
png(paste0(projectFolder, "/codex_XY_scatter_TLS_like_subset.png"), height=1200, width=1400)
print(FeatureScatter(Akoya_TLS_like_subset, 
               feature1 = "Centroid.X.µm", 
               feature2 = "Centroid.Y.µm", 
               raster = FALSE, 
               shuffle = TRUE, 
               group.by = "Name",
               pt.size = 6))
dev.off()

# Visualizing two scores blended
png(paste0(projectFolder, "imagefeatureplots_blended/codex_sample_P21_TLS_like_subset_blended.png"), height=1200, width=1400)
p1 <- ImageFeaturePlot(Akoya_TLS_like_subset, features = c("TLS_like_TME_score1","Type1_immunity_score1"), min.cutoff = "q05", max.cutoff = "q75", size = 6, alpha = 0.9, border.color = "NA", dark.background = F, blend = TRUE, combine = FALSE)
p1[[3]] + NoLegend()  # Get just the co-expression plot
dev.off()

png(paste0(projectFolder, "imagefeatureplots_blended/codex_sample_P21_TLS_like_subset_unique_blended.png"), height=1200, width=1400)
p1 <- ImageFeaturePlot(Akoya_TLS_like_subset, features = c("TLS_like_TME_markers_unique1","Type1_immunity_markers_unique1"), min.cutoff = "q05", max.cutoff = "q75", size = 6, alpha = 0.9, border.color = "NA", dark.background = F, blend = TRUE, combine = FALSE)
p1[[3]] + NoLegend()  # Get just the co-expression plot
dev.off()





###### Part 3.2: Focus on Active type1 example section ###### 

Akoya_Type1_immunity_subset <- subset(Akoya_Seurat_sample1, Centroid.X.µm > 2750 & Centroid.X.µm < 3410 & Centroid.Y.µm > 7050 & Centroid.Y.µm < 7600)

Akoya_Type1_immunity_subset
# An object of class Seurat
# 60 features across 5584 samples within 1 assay
# Active assay: Akoya_Mean (60 features, 60 variable features)
#  3 layers present: counts, data, scale.data
#  1 spatial field of view present: Akoya

# Plot XY scatter of all cells
png(paste0(projectFolder, "/codex_XY_scatter_P21_Type1_immunity_subset.png"), height=1200, width=1400)
print(FeatureScatter(Akoya_Type1_immunity_subset, 
               feature1 = "Centroid.X.µm", 
               feature2 = "Centroid.Y.µm", 
               raster = FALSE, 
               shuffle = TRUE, 
               group.by = "Name",
               pt.size = 6))
dev.off()

# Optional: reverse y-coordinates for visualization
Akoya_Type1_immunity_subset@meta.data$`Centroid.Y.µm` <- -Akoya_Type1_immunity_subset@meta.data$`Centroid.Y.µm` #to change meta.data 
Akoya_Type1_immunity_subset@images$Akoya@boundaries$centroids@coords[,1] <- -Akoya_Type1_immunity_subset@images$Akoya@boundaries$centroids@coords[,1] #to change FOV 

# Visualizing two scores blended
png(paste0(projectFolder, "imagefeatureplots_blended/codex_sample_P21_Type1_immunity_subset_blended.png"), height=1200, width=1400)
p1 <- ImageFeaturePlot(Akoya_Type1_immunity_subset, features = c("TLS_like_TME_score1","Type1_immunity_score1"), min.cutoff = "q05", max.cutoff = "q75", size = 6, alpha = 0.9, border.color = "NA", dark.background = F, blend = TRUE, combine = FALSE)
p1[[3]] + NoLegend()  # Get just the co-expression plot
dev.off()

png(paste0(projectFolder, "imagefeatureplots_blended/codex_sample_P21_Type1_immunity_subset_unique_blended.png"), height=1200, width=1400)
p1 <- ImageFeaturePlot(Akoya_Type1_immunity_subset, features = c("TLS_like_TME_markers_unique1","Type1_immunity_markers_unique1"), min.cutoff = "q05", max.cutoff = "q75", size = 6, alpha = 0.9, border.color = "NA", dark.background = F, blend = TRUE, combine = FALSE)
p1[[3]] + NoLegend()  # Get just the co-expression plot
dev.off()




