### Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(GSEABase)
library(tidyverse)

##==============================================================================
## Create a merged object from the T-cells of all cancer types
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

### Subset object to keep only T-cells (including the NK)
Idents(object) <- "Majorcelltype_annotation"
object <- subset(object, idents = "T-cell")
object



##====================================================================================
## Perform an initial analysis to identify and remove low quality and doublet clusters
##====================================================================================

#### Initialize Seurat object
object <- NormalizeData(objectNew, verbose = T)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = T) 
object <- ScaleData(object, verbose = T, vars.to.regress = c("orig.ident", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "TumorType", "Stress.Score1", "IFN.Score1"))
object <- RunPCA(object, features = VariableFeatures(object), npcs = 50, verbose = T)

### Determine statistical significant number of PCs with elbow plot
pdf("ElbowPlot.pdf")
ElbowPlot(object, ndims = 50)
dev.off()
# Number of PCs selected according to elbow plot: 18

### Run Harmony to correct for the technological batch effect
pdf("Harmony.pdf")
options(repr.plot.height = 2.5, repr.plot.width = 6)
object <- object %>% 
  RunHarmony(group.by.vars = "Technology", plot_convergence = TRUE)
dev.off()

### Compute the UMAP and determine clusters with Louvain algorithm
res = c(0.2, 0.4, 0.5, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8)
object <- object %>% 
  RunUMAP(reduction = "harmony", dims = 1:18) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:18) 
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

### Find differentially expressed genes among clusters according to chosen resolution (res.1) 
Idents(object) <- object$RNA_snn_res.1
DefaultAssay(object) <- "RNA"
object.markers <- FindAllMarkers(object, min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.25, return.thresh = 0.01, only.pos = FALSE)
write.table(object.markers, "object_markers_Res1.txt",col.names=NA, sep="\t")

### Plot marker genes to identify T-cell subtypes (see http://apps.lambrechtslab.sites.vib.be/PanCancer-Atlas/, under "Single-Cell Profiling" tab) 
DefaultAssay(object) <- "RNA"
Idents(object) <- "RNA_snn_res.1"
pdf("Markergenes_Tcells_subtypes.pdf")
genes <- c("MarkersToAdd")   
for(i in 1:length(genes)){
  print(FeaturePlot(object, features = genes[i], cols = c("grey", "blue")))
}
dev.off()

# Check if marker genes that identify other major cell types are expressed 
# - T-cell: "CD3D", "CD3E", "CD4", "CD8A"
# - B-cell: "CD79A"
# - Macrophages/monocytes: "CD68", "LYZ", "AIF1", "ITGAM"
# - pDC: "LILRA4", "CXCR3"
# - cDC: "CLEC9A", "XCR1", "CD1C", "CCR7", "CCL17", "CCL19"
# - Mast cell: "MS4A2", "TPSAB1", "CPA3"
# - Endothelial cell: "CLDN5", "PECAM1", "VWF"
# - Fibroblast: "COL1A1", "BGN", "DCN"
# - Epithelial cell: "EPCAM", "KRT7", "KRT18"

DefaultAssay(object) <- "RNA"
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

### Cluster Minor_celltype_annotation 
object <- RenameIdents(object, 
                  "0" = "T-cell",
                  "1" = "T-cell",
                  "2" = "T-cell",
                  "3" = "T-cell",
                  "4" = "T-cell",
                  "5" = "T-cell",
                  "6" = "T-cell",
                  "7" = "T-cell",
                  "8" = "Low quality", 
                  "9" = "NK cell",
                  "10" = "IFN",
                  "11" = "T-cell",
                  "12" = "NK cell",
                  "13" = "T-cell",
                  "14" = "Proliferative",
                  "15" = "T-cell",
                  "16" = "Doublets",
                  "17" = "T-cell",
                  "18" = "Low quality", 
                  "19" = "Doublets", 
                  "20" = "T-cell",
                  "21" = "T-cell",
                  "22" = "T-cell",
                  "23" = "T-cell",
                  "24" = "Doublets" 
)
object$RNA_snn_res.1 <- Idents(object)
object$Tcells_Minor_celltype_annotation <- object$RNA_snn_res.1

# We identified clusters of low quality cells (based of low nFeature and nCount, high percent.mt and very low gene expression) 
# and doublets (based on DoubletsFinder and expression of genes from multiple major cell types)

### Keep only the high-quality T or NK cells
object <- subset(object, idents=c("T-cell", "IFN", "NK cell", "Proliferative"))

### Minor_celltype_annotation visualization 
pdf("Minor_celltype_annotation.pdf")
Idents(object) <- "Tcells_Minor_celltype_annotation"
UMAPPlot(object)
dev.off()


##===========================================================================
## With the cleaned/high-quality object, repeat the above-mentioned analysis 
##===========================================================================

### Run the standard preprocessing
object <- NormalizeData(object)
object <- FindVariableFeatures(object, nfeatures=3000)

### Remove some variable features that may cause contamination in the data
vf <- VariableFeatures(object)[!grepl(pattern="^HB[AB]", VariableFeatures(object))]
VariableFeatures(object) <- vf[1:2000] # Keep top2000 non-contaminant genes

### Compute 50 PCs 
object <- ScaleData(object, vars.to.regress=regression_conditions)
object <- RunPCA(object, npcs=50)
ElbowPlot(object)

### We pick the top 40 PCs to continue with the analysis
object <- RunHarmony(object, 
                group.by.vars = "Technology", 
                plot_convergence = TRUE, 
                assay.use="RNA", 
                dims.use=1:40,
                epsilon.cluster = -Inf, # Force the harmony algorithm to not stop early
                epsilon.harmony = -Inf
)

### Perform clustering and 2D dimensionality reduction
resolutions <- c(seq(0, 2, 0.1), seq(2, 3, 0.2), seq(3.5, 6, 0.5))
object <- RunUMAP(object, dims=1:40, reduction="harmony")
object <- FindNeighbors(object, dims=1:40, reduction="harmony")
object <- FindClusters(object, resolution=resolutions)


### Manually annotate the clusters by looking at canonical gene markers (see http://apps.lambrechtslab.sites.vib.be/PanCancer-Atlas/, under "Single-Cell Profiling" tab) 
Idents(object) <- object$RNA_snn_res.0.7
object <- RenameIdents(object,
                  "0" = "CD4+ TEM",
                  "1" = "CD8+ TEM",
                  "2" = "CD4+ TREG",
                  "3" = "CD4+ TN",
                  "4" = "CD8+ TEX",
                  "5" = "CD8+ TRM",
                  "6" = "CD4+ TEX",
                  "7" = "NK cyto",
                  "8" = "NK inflam",
                  "9" = "CD4+ TEM",
                  "10" = "CD4+ TH17",
                  "11" = "Prolif T",
                  "12" = "CD8+ EMRA",
                  "13" = "CD8- γδ",
                  "14" = "Low quality",
                  "15" = "CD8+ TEM",
                  "16" = "Prolif T",
                  "17" = "CD8+ TEX",
                  "18" = "Low quality",
                  "19" = "MAIT",
                  "20" = "Low quality",
                  "21" = "Prolif T",
                  "22" = "Low quality"
)
object$Minor_celltype_annotation <- Idents(object)

# Additionally, we identify some rare Tcell subtypes which can only be separated at higher clustering resolutions
object$Minor_celltype_annotation[object$RNA_snn_res.1.8 == "18"] <- "CD8+ γδ"
object$Minor_celltype_annotation[object$RNA_snn_res.4.5 == "32"] <- "CD8+ TN"
object$Minor_celltype_annotation[object$RNA_snn_res.3.6 == "37"] <- "CD8+ MAIT"
object$annotation[object$RNA_snn_res.0.7 == "6"] <- "CD4+ TH1"
object$annotation[object$RNA_snn_res.2.2 == "13"] <- "CD4+ TFH"

### Remove low quality cluster
Ident(object) <- "Minor_celltype_annotation"
object <- subset(object, ident = "Low Quality", invert = T)
object

### Annotation visualization 
pdf("Annotation.pdf")
Idents(object) <- s$Minor_celltype_annotation
UMAPPlot(object)
dev.off()

### Save final Seurat object with annotation
saveRDS(object, file = "Tcells.rds")

### Due to a technical issue leading to sample misclassification, one HNSCC patient was excluded from further analysis to ensure the accuracy and reliability of our findings
object = subset(object, subset = Patient != "HNSCC_P13")
saveRDS(object, file = "Tcells.rds")