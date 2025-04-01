### Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(GSEABase)
library(tidyverse)


##==============================================================================
## Create a merged object from the B-cells of all cancer types
##==============================================================================

### Read object from individual cancer types (see file "Individual cancer types.R") 
BC_early <- readRDS("/path/BC_early.rds")
BC_advanced <- readRDS("/path/BC_advanced.rds")
CC <- readRDS("/path/CC.rds")
CRC <- readRDS("/path/CRC.rds")
GBM <- readRDS("/path/GBM.rds")
HCC <- readRDS("/path/HCC.rds")
HGSOC <- readRDS("/path/HGSOC.rds")
HNSCC <- readRDS("/path/HNSCC.rds")
MEL <- readRDS("/path/MEL.rds")
NSCLC_early <- readRDS("/path/NSCLC_early.rds")
NSCLC_advanced <- readRDS("/path/NSCLC_A.rds")

### Merge objects
object <- merge(BC_early, y=c(BC_advanced, CC, CRC, GBM, HCC, HGSOC, HNSCC, MEL, NSCLC_early, NSCLC_advanced))

### Subset object to keep only B-cells
Idents(object) <- "Majorcelltype_annotation"
object <- subset(object, idents = "B-cell")
object



##====================================================================================
## Perform an initial analysis to identify and remove low quality and doublet clusters
##====================================================================================

# Exclude contaminant features like B-cell receptor genes and Hemoglobin
IG.gene1 <- grep(pattern = "^IGLV", x = rownames(object), value = TRUE)
IG.gene2 <- grep(pattern = "^IGKV", x = rownames(object), value = TRUE)
IG.gene3 <- grep(pattern = "^IGHV", x = rownames(object), value = TRUE)
IG.genes <- c(IG.gene1, IG.gene2, IG.gene3)
hgGenes = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
excluded_features = c(IG.genes, hgGenes)

#### Initialize Seurat object
object <- object %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000 + length(excluded_features))
VariableFeatures(object) <- VariableFeatures(object)[!(VariableFeatures(s) %in% excluded_features)][1:2000]

object <- object %>%
  ScaleData(vars.to.regress = c("orig.ident", "TumorType3", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score")) %>%
  RunPCA(npcs=50) 

### Determine statistical significant number of PCs with elbow plot
pdf("ElbowPlot.pdf")
ElbowPlot(object, ndims = 50)
dev.off()
# Number of PCs selected according to elbow plot: 40

### Compute the UMAP and determine clusters with Louvain algorithm
pcs <- 40
object = object %>%
  RunUMAP(dims=1:pcs) %>%
  FindNeighbors(dims=1:pcs) %>% 
  FindClusters(resolution=seq(0, 1, 0.1))

### Plots
pdf("Resolutions.pdf")
for (i in 1:14){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(object, reduction = "umap", group.by = resolution) + ggtitle(resolution))
  print(DimPlot(object, reduction = "umap", group.by = resolution, label = TRUE) + ggtitle(resolution))
}
dev.off()

### DoubletFinder visualization
pdf("Doublets.pdf")
DimPlot(object, group.by = "pANNPredictions")
DimPlot(object, group.by = "pANNPredictions", reduction = "umap")
dev.off()

### Remove doublet cluster 
object <- subset(object, subset = RNA_snn_res.0.1 != "2")
object


##==============================================================================================================================
## Integrated the cleaned data with Seurat CCA; Harmony was found to not be strong enough to remove the batch effects in B-cells
##==============================================================================================================================


### Split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(object, split.by = "Technology")
nvf = 2000

### Normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = nvf + length(excluded_features))
  VariableFeatures(x) <- VariableFeatures(x)[!(VariableFeatures(x) %in% excluded_features)][1:nvf]
  return(x)
})

### Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

### Find anchors
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

### Dimensionality reduction and clustering
pcs = 12
combined <- ScaleData(combined, verbose = T, vars.to.regress = c("orig.ident", "nCount_RNA", "TumorType3", "percent.mt", "S.Score", "G2M.Score", "IFN.Score1", "Stress.Score1", "Hypoxia.Score1"))
combined <- RunPCA(combined, npcs = pcs, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:pcs)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:pcs)
combined <- FindClusters(combined, resolution = seq(0, 2, 0.1))

### Annotate B-cell clusters based on canonical gene markers (see http://apps.lambrechtslab.sites.vib.be/PanCancer-Atlas/, under "Single-Cell Profiling" tab) 
combined$Minor_celltype_annotation <- NA
combined@meta.data[combined$integrated_snn_res.0.8 %in% c(13, 15), "Minor_celltype_annotation"] <- "GC B"
combined@meta.data[combined$integrated_snn_res.0.8 %in% c(4, 10), "Minor_celltype_annotation"] <- "IgA mature"
combined@meta.data[combined$integrated_snn_res.0.8 %in% c(5, 7, 8), "Minor_celltype_annotation"] <- "IgG mature"
combined@meta.data[combined$integrated_snn_res.0.8 %in% c(6), "Minor_celltype_annotation"] <- "IgG immature"
combined@meta.data[combined$integrated_snn_res.0.8 %in% c(12), "Minor_celltype_annotation"] <- "IgG mature"
combined@meta.data[combined$integrated_snn_res.1 %in% c(2), "Minor_celltype_annotation"] <- "Naive mature"
combined@meta.data[combined$integrated_snn_res.1 %in% c(15), "Minor_celltype_annotation"] <- "Plasmablast"
combined@meta.data[combined$integrated_snn_res.1.9 %in% c(15), "Minor_celltype_annotation"] <- "IgA immature"
combined@meta.data[combined$integrated_snn_res.2 %in% c(3), "Minor_celltype_annotation"] <- "Breg"
combined@meta.data[combined$integrated_snn_res.2 %in% c(17, 24), "Minor_celltype_annotation"] <- "Memory IgM-"
combined@meta.data[is.na(combined$Minor_celltype_annotation), "Minor_celltype_annotation"] <- "Low quality"

pdf(file="CCA_Plots.pdf")
DimPlot(object, reduction = "umap", group.by = "orig.ident")
DimPlot(object, reduction = "umap", group.by = "Patient")
DimPlot(object, reduction = "umap", group.by = "Technology")
DimPlot(object, reduction = "umap", group.by = "TumorType")
FeaturePlot(object, features = "nFeature_RNA")
FeaturePlot(object, features = "nCount_RNA")
FeaturePlot(object, features = "percent.mt")
FeaturePlot(object, features = "S.Score")
FeaturePlot(object, features = "G2M.Score")
FeaturePlot(object, features = "MKI67")
FeaturePlot(object, features = "EMT.Score1")
FeaturePlot(object, features = "IFN.Score1")
FeaturePlot(object, features = "Stress.Score1")
FeaturePlot(object, features = "Hypoxia.Score1")
FeaturePlot(object, features = "IG.Score1")
dev.off()

### Plots
pdf("CCA_Resolutions.pdf")
for (i in 1:14){
  resolution <- paste0("integrated_snn_res." , res[i])
  print(DimPlot(object, reduction = "umap", group.by = resolution) + ggtitle(resolution))
  print(DimPlot(object, reduction = "umap", group.by = resolution, label = TRUE) + ggtitle(resolution))
}
dev.off()

# We identified clusters of low quality cells (based of low nFeature and nCount, high percent.mt and very low gene expression) 

### Remove low quality cluster
Ident(object) <- "Minor_celltype_annotation"
object <- subset(object, ident = "Low Quality", invert = T)
object

### Annotation visualization 
pdf("Annotation.pdf")
Idents(object) <- "Minor_celltype_annotation"
UMAPPlot(object)
dev.off()

### Save Seurat object with annotation
saveRDS(object, file = "Bcells.rds")

### Due to a technical issue leading to sample misclassification, one HNSCC patient was excluded from further analysis to ensure the accuracy and reliability of our findings
object = subset(object, subset = Patient != "HNSCC_P13")
saveRDS(object, file = "Bcells.rds")


