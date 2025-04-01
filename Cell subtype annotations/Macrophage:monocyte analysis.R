### Perform myeloid cell subclustering analysis across different cancer types with the sole aim of identifying cancer type-specific myeloid populations. 
# This approach enabled us to identify Kupffer cells in HCC, microglia in GBM and alveolar macrophages in NSCLC.
# The analysis followed the same steps as the one found in "Individual cancer types.R", using the following marker genes, collected from literature studies: 
# 1) Kupffer cells in HCC (MacParland et al, Nat. Commun. 2018; Aizarani et al, Nature 2019; Cappuyins et al, Nat. Commun. 2023): "CD163", "MAFB", "VSIG4", "LYVE1", "CD209", "PLTP", "SELENOP", "MERTK"
# 2) Microglia in GBM (Darmanis et al, Cell Rep. 2017; Masuda et al, Nature 2019; Muller et al, Genome Biology 2017): "TMEM119", "P2RY12", "P2RY13", "SLC2A5", "CX3CR1", "FCGR3B"
# 3) Alveolar macrophages in NSCLC (Qian et al, Cell Research 2020; Chen et al, Cell 2021): "FABP4", "MARCO", "PPARG", "INHBA", "ALDH2"

### Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(GSEABase)
library(tidyverse)

##==============================================================================
## Create a merged object from the macrophages/monocytes of all cancer types
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

### Subset objet to keep only macrophages/monocytes 
Idents(object) <- "Majorcelltype_annotation"
object <- subset(object, idents = "Macrophages/monocytes")
object
#76436 cells 

### Remove the cancer-specific macrophages/monocytes mentioned above 
# Kupffer cells (HCC)
kupffer = rownames(dplyr::filter(HCC_Myeloid@meta.data, MyeloidAnnotation == "Kupffer_cell"))
kupffer = intersect(kupffer, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% kupffer]) # Only select cells that are not kupffer cells

# Microglia (GBM)
microglia = rownames(dplyr::filter(GBM_Myeloid@meta.data, MyeloidAnnotation == "Microglia_cell"))
microglia = intersect(microglia, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% microglia]) # Only select cells that are not microglia cells

# Alveolar macrophages (NSCLC)
alveolar = rownames(dplyr::filter(NSCLC_Myeloid@meta.data, MyeloidAnnotation == "Alveolar_macrophages"))
alveolar = intersect(alveolar, colnames(object)) # Make sure that only common cells are selected
object = subset(object, cells = colnames(object)[!object %in% alveolar]) # Only select cells that are not alveolar macrophages

object 
# 60506 cells



##====================================================================================
## Perform an initial analysis to identify and remove low quality and doublet clusters
##====================================================================================

#### Initialize Seurat object
object <- NormalizeData(objectNew, verbose = T)
object <-FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = T) 
object <-ScaleData(object, verbose = T, vars.to.regress = c("orig.ident", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "TumorType", "Hypoxia.Score1", "IFN.Score1"))
object <-RunPCA(object, features = VariableFeatures(object), npcs = 50, verbose = T)

### Determine statistical significant number of PCs with elbow plot
pdf("ElbowPlot.pdf")
ElbowPlot(object, ndims = 50)
dev.off()
# Number of PCs selected according to elbow plot: 17 

### Run Harmony to correct for the technological batch effect
pdf("Harmony.pdf")
options(repr.plot.height = 2.5, repr.plot.width = 6)
object <- object %>% 
  RunHarmony("Technology", plot_convergence = TRUE)
dev.off()

### Compute the UMAP and determine clusters with Louvain algorithm
object <- object %>% 
  RunUMAP(reduction = "harmony", dims = 1:17) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:17) 
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

### Find differentially expressed genes among clusters according to chosen resolution (0.5) 
Idents(object) <- "RNA_snn_res.0.5"
DefaultAssay(object) <- "RNA"
object.markers_Res0.5 <- FindAllMarkers(object, min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.25, return.thresh = 0.01, only.pos = FALSE)
write.table(object.markers, "object_markers_Res0.5.txt",col.names=NA, sep="\t")

### Plot marker genes to identify Macrophages/monocytes subtypes (see http://apps.lambrechtslab.sites.vib.be/PanCancer-Atlas/, under "Single-Cell Profiling" tab) 
DefaultAssay(object) <- "RNA"
Idents(object) <- "RNA_snn_res.0.5"
pdf("Markergenes_Macrophages/monocytes_subtypes.pdf")
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
Idents(object) <- "RNA_snn_res.0.5"
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
Idents(object) <- "RNA_snn_res.0.5"
object <- RenameIdents(object,
                       `0` = "LAM1", `1` = "Low Quality", `2` = "Classical Mono",`3` = "Suppr Mac", `4` = "Hypoxic Mac", `5` = "IFN Mac",
                       `6` = "Hypoxic Mac", `7` = "LAM2", `8` = "MT Mac",  `9` = "Non-Classical Mono", `10` = "Low Quality",
                       `11` = "Neutrophils", `12` = "Prolif Mac", `13` = "Low Quality")

object$RNA_snn_res.0.5 <- Idents(object)
object$Minor_celltype_annotation <- object$RNA_snn_res.0.5

# We identified clusters of low quality cells (based of low nFeature and nCount, high percent.mt and very low gene expression) 

#### Add annotation of specific clusters at different resolution
object$Minor_celltype_annotation <- as.character(object$Minor_celltype_annotation)
Idents(object) <- "RNA_snn_res.1"
cells <- WhichCells(object, expression = RNA_snn_res.1 == "3")
object$Minor_celltype_annotation[cells] <- "Mono-like Mac"

Idents(object) <- "RNA_snn_res.0.8"
cells <- WhichCells(object, expression = RNA_snn_res.0.8 == "15")
object$Minor_celltype_annotation[cells] <- "Inflam Mac"
cells <- WhichCells(object, expression = RNA_snn_res.0.8 == "4")
object$Minor_celltype_annotation[cells] <- "Mono-like LAM"
cells <- WhichCells(object, expression = RNA_snn_res.0.8 == "3")
object$Minor_celltype_annotation[cells] <- "Perivasc Mac"
cells <- WhichCells(object, expression = RNA_snn_res.0.8 == "20")
object$Minor_celltype_annotation[cells] <- "Perivasc Mac"
cells <- WhichCells(object, expression = RNA_snn_res.0.8 == "6")
object$Minor_celltype_annotation[cells] <- "Perivasc Mac"

Idents(object) <- "RNA_snn_res.1.6"
cells <- WhichCells(object, expression = RNA_snn_res.1.6 == "19")
object$Minor_celltype_annotation[cells] <- "Interm Mac"

object$Minor_celltype_annotation <- as.factor(object$Minor_celltype_annotation)

### Remove low quality cluster
Ident(object) <- "Minor_celltype_annotation"
object <- subset(object, ident = "Low Quality", invert = T)
object

### Annotation visualization 
pdf("Annotation.pdf")
Idents(object) <- "Minor_celltype_annotation"
UMAPPlot(object)
dev.off()

object
# 46150 cells

### save Seurat object with annotation
saveRDS(object, file = "MacrophagesMonocytes.rds")

### Due to a technical issue leading to sample misclassification, one HNSCC patient was excluded from further analysis to ensure the accuracy and reliability of our findings
object = subset(object, subset = Patient != "HNSCC_P13")
saveRDS(object, file = "MacrophagesMonocytes.rds")
