### Perform endothelial cell (EC) subclustering analysis across individual cancer types with the sole aim of identifying cancer type-specific EC populations. 
# This approach enabled us to identify alveolar ECs in NSCLC and liver sinusoidal ECs in HCC. 
# The analysis followed the same steps as the one found in "Individual cancer types.R", using the following marker genes, collected from literature studies: 
# 1) Alveolar ECs in NSCLC (Qian et al, Cell Research 2020; Astrid et al, Nature 2020): "EDNRB","IL1RL1","HPGD", "CLDN18", "AGER", "SFTPC"
# 2) Liver sinusoidal ECs in HCC (Aizarani et al, Nature 2019; MacParland et al, Nat. Commun. 2018): "CLEC4G","CLEC4M","FLT1","PECAM1","VWF","CD34", "MGP", "CCL14",	"DNASE1L3", "SPARCL1",	"CLEC1B",	"LIFR", "TM4SF1",	"FCN2",	"PTGDS", "CLEC14A",	"S100A13",	"C7"

### Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(GSEABase)
library(tidyverse)

##==============================================================================
## Create a merged object from the ECs of all cancer types
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

### Subset objet to keep only ECs
Idents(object) <- "Majorcelltype_annotation"
object <- subset(object, idents = "EC")
object
#26896 cells 


### Remove the cancer-specific ECs mentioned above 
# Alveolar ECs (NSCLC)
Alveolar = rownames(dplyr::filter(NSCLC_ECs@meta.data, ECAnnotation == "Alveolar_cell"))
Alveolar = intersect(Alveolar, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% Alveolar]) # Only select cells that are not Alveolar cells

# Liver sinusoidal ECs (HCC)
Sinusoidal = rownames(dplyr::filter(HCC_ECs@meta.data, ECAnnotation == "Sinusoidal_cell"))
Sinusoidal = intersect(Sinusoidal, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% Sinusoidal]) # Only select cells that are not Sinusoidal cells

# LQ and doublets 
LQ_NSCLC_EC = rownames(dplyr::filter(NSCLC_ECs@meta.data, ECAnnotation == "Low quality"))
LQ_NSCLC_EC = intersect(LQ_NSCLC_EC, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% LQ_NSCLC_EC]) # Only select cells that are not LQ_NSCLC_EC cells

Doublets_NSCLC_EC = rownames(dplyr::filter(NSCLC_ECs@meta.data, ECAnnotation == "Doublets"))
Doublets_NSCLC_EC = intersect(Doublets_NSCLC_EC, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% Doublets_NSCLC_EC]) # Only select cells that are not Doublets_NSCLC_EC cells

LQ_HCC_EC = rownames(dplyr::filter(HCC_ECs@meta.data, ECAnnotation == "Low quality"))
LQ_HCC_EC = intersect(LQ_HCC_EC, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% LQ_HCC_EC]) # Only select cells that are not LQ_HCC_EC cells

Doublets_HCC_EC = rownames(dplyr::filter(HCC_ECs@meta.data, ECAnnotation == "Doublets"))
Doublets_HCC_EC = intersect(Doublets_HCC_EC, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% Doublets_HCC_EC]) # Only select cells that are not Doublets_HCC_EC cells

object 
# 23509 cells



##====================================================================================
## Perform an initial analysis to identify and remove low quality and doublet clusters
##====================================================================================

#### Initialize Seurat object
object <- NormalizeData(objectNew, verbose = T)
object <-FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = T) 
object <-ScaleData(object, verbose = T, vars.to.regress = c("orig.ident", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "TumorType", "Stress.Score1", "IFN.Score1"))
object <-RunPCA(object, features = VariableFeatures(object), npcs = 50, verbose = T)

### Determine statistical significant number of PCs with elbow plot
pdf("ElbowPlot.pdf")
ElbowPlot(object, ndims = 50)
dev.off()
# Number of PCs selected according to elbow plot: 20 

### Run Harmony to correct for the technological batch effect
pdf("Harmony.pdf")
options(repr.plot.height = 2.5, repr.plot.width = 6)
object <- object %>% 
  RunHarmony("Technology", plot_convergence = TRUE)
dev.off()

### Compute the UMAP and determine clusters with Louvain algorithm
object <- object %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) 
res = c(0.2, 0.4, 0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8)
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

# Find differentially expressed genes among clusters according to chosen resolution (0.4) 
Idents(object) <- "RNA_snn_res.0.4"
DefaultAssay(object) <- "RNA"
object.markers_Res0.4 <- FindAllMarkers(object, min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.25, return.thresh = 0.01, only.pos = FALSE)
write.table(object.markers, "object_markers_Res0.4.txt",col.names=NA, sep="\t")

### Plot marker genes to identify EC subtypes (see http://apps.lambrechtslab.sites.vib.be/PanCancer-Atlas/, under "Single-Cell Profiling" tab) 
DefaultAssay(object) <- "RNA"
Idents(object) <- "RNA_snn_res.0.4"
pdf("Markergenes_EC_subtypes.pdf")
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
Idents(object) <- "RNA_snn_res.0.4"
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
Idents(object) <- "RNA_snn_res.0.4"
object <- RenameIdents(object,
                       `0` = "Low Quality", `1` = "PCV", `2` = "Stalk",`3` = "Arterial", `4` = "Tip", `5` = "Lymphatic",
                       `6` = "Capillary", `7` = "IFN EC", `8` = "Doublets",  `9` = "Doublets", `10` = "Doublets",`11` = "Doublets", 
                       `12` = "Doublets", `13` = "Prolif EC")
object$RNA_snn_res.0.4 <- Idents(object)
object$Minor_celltype_annotation <- object$RNA_snn_res.0.4

# We identified clusters of low quality cells (based of low nFeature and nCount, high percent.mt and very low gene expression) 
# and doublets (based on DoubletsFinder and expression of genes from multiple major cell types)

# Annotation visualization 
pdf("Annotation.pdf")
Idents(object) <- "Minor_celltype_annotation"
UMAPPlot(object)
dev.off()

### Remove low quality and doublet clusters
Idents(object) <- "Minor_celltype_annotation"
object <- subset(object, features = c("Doublets", "Low Quality"), inv = T)

object
# 15910 cells


##===========================================================================
## With the cleaned/high-quality object, repeat the above-mentioned analysis 
##===========================================================================

object <- NormalizeData(object, verbose = T)
object <-FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = T) 
object <-ScaleData(object, verbose = T, vars.to.regress = c("orig.ident", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "TumorType", "Stress.Score1", "IFN.Score1"))
object <-RunPCA(object, features = VariableFeatures(object), npcs = 50, verbose = T)

### Determine statistical significant number of PCs with elbow plot
pdf("ElbowPlot.pdf")
ElbowPlot(object, ndims = 50)
dev.off()
# Number of PCs selected according to elbow plot: 13 

### Run Harmony to correct for the technological batch effect
pdf("Harmony.pdf")
options(repr.plot.height = 2.5, repr.plot.width = 6)
object <- object %>% 
  RunHarmony("Technology", plot_convergence = TRUE)
dev.off()

### Compute the UMAP and determine clusters with Louvain algorithm
object <- object %>% 
  RunUMAP(reduction = "harmony", dims = 1:13) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:13) 
res = c(0.2, 0.4, 0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8)
object <- FindClusters(object, resolution = res) 

### Cluster annotation 
Idents(object) <- "RNA_snn_res.0.4"
object <- RenameIdents(object,
                       `0` = "Low Quality", `1` = "Stalk", `2` = "PCV",`3` = "Arterial", `4` = "Tip", `5` = "Lymphatic",
                       `6` = "Capillary", `7` = "IFN EC", `8` = "Prolif EC")
object$RNA_snn_res.0.4 <- Idents(object)
object$Minor_celltype_annotation <- object$RNA_snn_res.0.4

### Add annotation of specific clusters at different resolution to distinguish PCV and venous ECs 
object$Minor_celltype_annotation <- as.character(object$Minor_celltype_annotation)
cells <- WhichCells(object, expression = RNA_snn_res.1.6 == "5")
object$Minor_celltype_annotation[cells] <- "Venous"
object$Minor_celltype_annotation <- as.factor(object$Minor_celltype_annotation)

# We identified clusters of low quality cells (based of low nFeature and nCount, high percent.mt and very low gene expression) 

### Remove low quality cluster
Ident(object) <- "Minor_celltype_annotation"
object <- subset(object, ident = "Low Quality", invert = T)
object

### Annotation visualization 
pdf("Annotation_Highquality.pdf")
Idents(object) <- "Minor_celltype_annotation"
UMAPPlot(object)
dev.off()

### save Seurat object with annotation
saveRDS(object, file = "EC.rds")

### Due to a technical issue leading to sample misclassification, one HNSCC patient was excluded from further analysis to ensure the accuracy and reliability of our findings
object = subset(object, subset = Patient != "HNSCC_P13")
saveRDS(object, file = "DC.rds")