### Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(GSEABase)
library(tidyverse)

##==============================================================================
## Create a merged object from the DCs of all cancer types
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
NSCLC_advanced <- readRDS("/path/.rds")

### Merge objects
object <- merge(BC_early, y=c(BC_advanced, CC, CRC, GBM, HCC, HGSOC, HNSCC, MEL, NSCLC_early, NSCLC_advanced))

### Subset objet to keep only dendritic cells (DCs) 
Idents(object) <- "Majorcelltype_annotation"
object <- subset(object, idents = c("pDC", "cDC"))
object
#9806 cells 


##====================================================================================
## Perform an initial analysis to identify and remove low quality and doublet clusters
##====================================================================================

#### Initialize Seurat object
object <- NormalizeData(objectNew, verbose = T)
object <-FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = T) 
object <-ScaleData(object, verbose = T, vars.to.regress = c("orig.ident", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "TumorType"))
object <-RunPCA(object, features = VariableFeatures(object), npcs = 50, verbose = T)

### Determine statistical significant number of PCs with elbow plot
pdf("ElbowPlot.pdf")
ElbowPlot(object, ndims = 50)
dev.off()
# Number of PCs selected according to elbow plot: 11

### Run Harmony to correct for the technological batch effect
pdf("Harmony.pdf")
options(repr.plot.height = 2.5, repr.plot.width = 6)
object <- object %>% 
  RunHarmony("Technology", plot_convergence = TRUE)
dev.off()

### Compute the UMAP and determine clusters with Louvain algorithm
object <- object %>% 
  RunUMAP(reduction = "harmony", dims = 1:11) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:11) 
res = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8)
object <- FindClusters(object, resolution = res) 

### Plots
pdf("Resolutions.pdf")
for (i in 1:14){
  resolution <- paste0("RNA_snn_res." , res[i])
  print(DimPlot(object, reduction = "umap", group.by = resolution) + ggtitle(resolution))
  print(DimPlot(object, reduction = "umap", group.by = resolution, label = TRUE) + ggtitle(resolution))
}
dev.off()

pdf(file="Plots.pdf")
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

### Find differentially expressed genes among clusters according to chosen resolution (0.6) 
Idents(object) <- "RNA_snn_res.0.6"
DefaultAssay(object) <- "RNA"
object.markers_Res0.6 <- FindAllMarkers(object, min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.25, return.thresh = 0.01, only.pos = FALSE)
write.table(object.markers, "object_markers_Res.0.6.txt",col.names=NA, sep="\t")

### Plot marker genes to identify DC subtypes (see http://apps.lambrechtslab.sites.vib.be/PanCancer-Atlas/, under "Single-Cell Profiling" tab) 
DefaultAssay(object) <- "RNA"
Idents(object) <- "RNA_snn_res.0.6"
pdf("Markergenes_DC_subtypes.pdf")
genes <- c("MarkersToAdd")   
for(i in 1:length(genes)){
  print(FeaturePlot(object, features = genes[i], cols = c("grey", "blue")))
}
dev.off()

# Check if marker genes that identify other major cell types are expressed 
# T-cell: "CD3D", "CD3E", "CD4", "CD8A"
# B-cell: "CD79A"
# Macrophages/monocytes: "CD68", "LYZ", "AIF1", "ITGAM"
# pDC: "LILRA4", "CXCR3"
# cDC: "CLEC9A", "XCR1", "CD1C", "CCR7", "CCL17", "CCL19"
# Mast cell: "MS4A2", "TPSAB1", "CPA3"
# Endothelial cell: "CLDN5", "PECAM1", "VWF"
# Fibroblast: "COL1A1", "BGN", "DCN"
# Epithelial cell: "EPCAM", "KRT7", "KRT18"

DefaultAssay(object) <- "RNA"
Idents(object) <- "RNA_snn_res.0.6"
pdf("Markergenes_Majorcelltypes.pdf")
genes <- c("CD3D", "CD3E", "CD4", "CD8A", "CD79A", "CD68", "LYZ", "AIF1", "ITGAM", "LILRA4", "CXCR3", "CLEC9A", "XCR1", "CD1C", "CCR7", "CCL17", "CCL19",
           "MS4A2", "TPSAB1", "CPA3", "CLDN5", "PECAM1", "VWF", "COL1A1", "BGN", "DCN", "EPCAM", "KRT7", "KRT18")   
for(i in 1:length(genes)){
  print(FeaturePlot(object, features = genes[i], cols = c("grey", "blue")))
}
dev.off()

### DoubletFinder visualization 
pdf("Doublets.pdf")
DimPlot(object, group.by = "pANNPredictions")
DimPlot(object, group.by = "pANNPredictions", reduction = "umap")
dev.off()

### Cluster annotation 
Idents(object) <- "RNA_snn_res.0.6"
object <- RenameIdents(object,
                       "0" = "cDC2", "1" = "cDC2", "2" = "pDC", "3" = "pDC", "4" = "mQuiescDC", "5" = "mRegDC", "6" = "Lang-likeDC",
                       "7" = "cDC1", "8" = "cDC2", "9" = "Doublets", "10" = "pDC", "11" = "pDC", "12" = "Doublets", "13" = "pDC")
object$RNA_snn_res.0.6 <- Idents(object)
object$Minor_celltype_annotation <- object$RNA_snn_res.0.6

# We identified clusters of doublets (based on DoubletsFinder and expression of genes from multiple major cell types)

### Remove low quality cluster
Ident(object) <- "Minor_celltype_annotation"
object <- subset(object, ident = "Doublets", invert = T)
object

# Annotation visualization 
pdf("Annotation.pdf")
Idents(object) <- "Minor_celltype_annotation"
UMAPPlot(object)
dev.off()

object
# 9460 cells

### save Seurat Object with annotation
saveRDS(object, file = "DC.rds")

### Due to a technical issue leading to sample misclassification, one HNSCC patient was excluded from further analysis to ensure the accuracy and reliability of our findings
object = subset(object, subset = Patient != "HNSCC_P13")
saveRDS(object, file = "DC.rds")
