# Load libraries 
library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(GSEABase)
library(tidyverse)
library(Matrix)
library(dplyr)
source("/Master_files/ReadEnsemble.R")
source("/Master_files/EnsembleToGenes.R")

# When available (i.e., for BC_early, CRC, HCC, HNSCC and NSCLC_early), Seurat objects with annotation of the main cell types from previous studies of our lab were used as starting point of the analysis:
# BC_early: Bassez et al, Nat. Med. 2021  
# CRC: Qian et al, Cell Research 2020
# HCC: Cappuyns et al, Nat. Commun 2023
# HNSCC: Franken et al, Immunity 2024
# NSCLC_early: Lambrecths et al, Nat. Med. 2018

# For all the other cancer types, we performed new analysis to identify major cell types. The codes below use HGSOC as representative example.

### CellRanger (to repeat for each sample)
# cellranger count --id= OV_1 \
# --transcriptome=/opt/refdata-gex-GRCh38-2020-A \
# --fastqs=/fastq_path \
# --sample=mysample \
# --create-bam=true \
# --localcores=8 \
# --localmem=64

### Setup the Seurat object
# Sample OV_1
OV_1 <-  ReadEnsemble("/OV_1/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_1), value=T, invert=T)

OV_1 <-  EnsembleToGenes (OV_1)
grep ("ENSG" , rownames(OV_1), value=T)

OV_1 <- CreateSeuratObject(counts = OV_1 ,project = "OV_1", min.cells = 1, min.features = 1)

# Sample OV_2
OV_2 <-  ReadEnsemble("/OV_2/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_2), value=T, invert=T)

OV_2 <-  EnsembleToGenes (OV_2)
grep ("ENSG" , rownames(OV_2), value=T)

OV_2 <- CreateSeuratObject(counts = OV_2 ,project = "OV_2", min.cells = 1, min.features = 1)

# Sample OV_3
OV_3 <-  ReadEnsemble("/OV_3/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_3), value=T, invert=T)

OV_3 <-  EnsembleToGenes (OV_3)
grep ("ENSG" , rownames(OV_3), value=T)

OV_3 <- CreateSeuratObject(counts = OV_3 ,project = "OV_3", min.cells = 1, min.features = 1)

# Sample OV_4
OV_4 <-  ReadEnsemble("/OV_4/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_4), value=T, invert=T)

OV_4 <-  EnsembleToGenes (OV_4)
grep ("ENSG" , rownames(OV_4), value=T)

OV_4 <- CreateSeuratObject(counts = OV_4 ,project = "OV_4", min.cells = 1, min.features = 1)

# Sample OV_5
OV_5 <-  ReadEnsemble("/OV_5/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_5), value=T, invert=T)

OV_5 <-  EnsembleToGenes (OV_5)
grep ("ENSG" , rownames(OV_5), value=T)

OV_5 <- CreateSeuratObject(counts = OV_5 ,project = "OV_5", min.cells = 1, min.features = 1)

# Sample OV_6
OV_6 <-  ReadEnsemble("/OV_6/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_6), value=T, invert=T)

OV_6 <-  EnsembleToGenes (OV_6)
grep ("ENSG" , rownames(OV_6), value=T)

OV_6 <- CreateSeuratObject(counts = OV_6 ,project = "OV_6", min.cells = 1, min.features = 1)

# Sample OV_7
OV_7 <-  ReadEnsemble("/OV_7/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_7), value=T, invert=T)

OV_7 <-  EnsembleToGenes (OV_7)
grep ("ENSG" , rownames(OV_7), value=T)

OV_7 <- CreateSeuratObject(counts = OV_7 ,project = "OV_7", min.cells = 1, min.features = 1)

# Sample OV_8
OV_8 <-  ReadEnsemble("/OV_8/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_8), value=T, invert=T)

OV_8 <-  EnsembleToGenes (OV_8)
grep ("ENSG" , rownames(OV_8), value=T)

OV_8 <- CreateSeuratObject(counts = OV_8 ,project = "OV_8", min.cells = 1, min.features = 1)

# Sample OV_9
OV_9 <-  ReadEnsemble("/OV_9/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_9), value=T, invert=T)

OV_9 <-  EnsembleToGenes (OV_9)
grep ("ENSG" , rownames(OV_9), value=T)

OV_9 <- CreateSeuratObject(counts = OV_9 ,project = "OV_9", min.cells = 1, min.features = 1)

# Sample OV_10
OV_10 <-  ReadEnsemble("/OV_10/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_10), value=T, invert=T)

OV_10 <-  EnsembleToGenes (OV_10)
grep ("ENSG" , rownames(OV_10), value=T)

OV_10 <- CreateSeuratObject(counts = OV_10 ,project = "OV_10", min.cells = 1, min.features = 1)

# Sample OV_11
OV_11 <-  ReadEnsemble("/OV_11/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_11), value=T, invert=T)

OV_11 <-  EnsembleToGenes (OV_11)
grep ("ENSG" , rownames(OV_11), value=T)

OV_11 <- CreateSeuratObject(counts = OV_11 ,project = "OV_11", min.cells = 1, min.features = 1)

# Sample OV_12
OV_12 <-  ReadEnsemble("/OV_12/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_12), value=T, invert=T)

OV_12 <-  EnsembleToGenes (OV_12)
grep ("ENSG" , rownames(OV_12), value=T)

OV_12 <- CreateSeuratObject(counts = OV_12 ,project = "OV_12", min.cells = 1, min.features = 1)

# Sample OV_13
OV_13 <-  ReadEnsemble("/OV_13/raw_gene_bc_matrices/GRCh38") 
grep ("ENSG" , rownames(OV_13), value=T, invert=T)

OV_13 <-  EnsembleToGenes (OV_13)
grep ("ENSG" , rownames(OV_13), value=T)

OV_13 <- CreateSeuratObject(counts = OV_13 ,project = "OV_13", min.cells = 1, min.features = 1)


### Merge samples
object <- merge(OV_1, y = c(OV_1, OV_2, OV_3, OV_4, OV_5, OV_6, OV_7, OV_8, OV_9, OV_10, OV_11, OV_12, OV_13), project = "HGSOC")

### Save Seurat object
# saveRDS(object, "/path/HGSOC.rds")

### Read Seurat object 
# object <- readRDS("/path/HGSOC.rds")

### Add percent.mt 
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

### Filtering n1: Default filtering, but without upper limit of nFeature (to keep the putative doublets)
object <- subset(object, subset = nFeature_RNA > 200 & percent.mt < 25 & nCount_RNA > 400 )

### Doublet finder
# orig.ident = individual samples
library(DoubletFinder)

Idents(object) <- "orig.ident"
object$pANN <- "NA"
object$pANNPredictions <- "NA"

for(sample in unique(object$orig.ident)){
  sample.cluster <- subset(object, idents = sample)
  print(paste0("sample:", sample))
  length(rownames(sample.cluster@meta.data))
  expected.doublets <- ceiling(0.039 * length(rownames(sample.cluster@meta.data)))
  sample.cluster <- doubletFinder_v3(sample.cluster, PCs = 1:20, nExp = expected.doublets, pN = 0.25, pK = 0.01)
  sample.cluster$pANN <- sample.cluster@meta.data[colnames(sample.cluster), paste("pANN_0.25_0.01", expected.doublets, sep = "_")]
  sample.cluster$pANNPredictions <- sample.cluster@meta.data[colnames(sample.cluster), paste("DF.classifications_0.25_0.01", expected.doublets, sep = "_")]
  object$pANN[colnames(sample.cluster)] <- sample.cluster$pANN[colnames(sample.cluster)]
  object$pANNPredictions[colnames(sample.cluster)] <- sample.cluster$pANNPredictions[colnames(sample.cluster)]
  sample.cluster <- NULL
}

### Filtering n2: Same as Filtering n1, but this time with upper limit of nFeature as well
# cancer types that required slightly different (more stringent in 3/4 cases) cutoff due to cancer-specific intrinsic properties: 
# GBM: nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20 & nCount_RNA > 600
# HCC: nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50 & nCount_RNA > 400
# HNSCC and BC early: nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15 & nCount_RNA > 400
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25 & nCount_RNA > 400 )

### Data normalization 
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose=T)

### Identification of highly variable features 
object <- FindVariableFeatures(object = object, selection.method = "vst", nfeatures = 2000)

### Scoring
# Files of specific gene signatures can be provived upon request
# Cell-cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
object <- CellCycleScoring(object = object, s.features = s.genes, g2m.features = g2m.genes)

# Interferon scoring
msigdb <- getBroadSets("/path/msigdb_v7.0.xml")
index <- sapply(msigdb, function(gs)
  bcCategory(collectionType(gs))=="c2")
geneset.c2 = msigdb[index]
geneset.interferon <- geneset.c2[["BROWNE_INTERFERON_RESPONSIVE_GENES"]]
interferon.genes <- geneset.interferon@geneIds[ geneset.interferon@geneIds %in% rownames(GetAssayData(object))]
object <- AddModuleScore(object = object, features = list(interferon.genes), 
                             ctrl = length(interferon.genes), name = 'IFN.Score')

# Hypoxia scoring
genes.hypoxia<- c("VEGFA", "SLC2A1", "PGAM1", "ENO1","LDHA", "TPI1", "P4HA1", "MRPS17",
                  "CDKN3", "ADM", "NDRG1", "TUBB6","ALDOA", "MIF", "ACOT7")
object <- AddModuleScore(object = object, features = list(genes.hypoxia),
                             ctrl = length(genes.hypoxia), name = 'Hypoxia.Score')

# IG scoring
IGH <- grep("^IGH", VariableFeatures(object), value = TRUE)
IGL <- grep("^IGL", VariableFeatures(object), value = TRUE)
IGK <- grep("^IGK", VariableFeatures(object), value = TRUE)
IGH <- grep("^IGH", rownames(object), value = TRUE)
IGL <- grep("^IGL", rownames(object), value = TRUE)
IGK <- grep("^IGK", rownames(object), value = TRUE)
IG_genes <- c(IGH, IGL, IGK)
object <- AddModuleScore(object = object, features = list(IG_genes), ctrl = length(IG_genes), name = 'IG.Score')

# Stress scoring
stress.genes <- read_csv("/path/stress.genes.csv") %>% .$gene
stress.genes <- stress.genes[stress.genes %in% rownames(GetAssayData(object))]
object <- AddModuleScore(object = object, features = list(stress.genes),
                             ctrl = length(stress.genes), name = 'Stress.Score')

# EMT scoring
msigdb.emt <- getBroadSets("/path/geneset_emt.xml")
index <- sapply(msigdb.emt, function(gs)
  bcCategory(collectionType(gs))=="c2")
geneset = msigdb.emt[index]
geneset.emt <- geneset[["SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN"]]
emt.genes <- geneset.emt@geneIds[ geneset.emt@geneIds %in% rownames(GetAssayData(object))]
object <- AddModuleScore(object = object, features = list(emt.genes), 
                             ctrl = length(emt.genes), name = 'EMT.Score')

# Epithelial differentiation scoring
genes.EpithelialDiff <- c("EPCAM", "KRT16", "KRT6A", "KRT6B", "KRT6C", "CLDN4", "KLK5", "KLK6", "KLK7", "KLK8", "KLK9", "KLK10", "KLK11", "KRT17", "KRT75", "S100A7", "S100A8", "S100A9", "CLDN1", "CLDN7", "SPRR1B")
object <- AddModuleScore(object = object, features = list(genes.EpithelialDiff),
                             ctrl = length(genes.EpithelialDiff), name = 'EpithelialDiff.Score')

### Data scaling
object <- ScaleData(object, verbose=T, vars.to.regress = c("orig.ident", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

### Determine statistical significant number of PCs with elbow plot
object <- RunPCA(object, features = VariableFeatures(object = object), npcs = 50, ndims.print = 1:5, nfeatures.print = 30)

pdf("ElbowPlot.pdf")
ElbowPlot(object, ndims = 50)
dev.off()

# Number of PCs selected according to elbow plot per cancer type: 
# BC early, nPC = 20
# BC advanced, nPC = 14
# CC, nPC = 27
# CRC, nPC = 30
# GBM, nPC = 14
# HCC, nPC = 12
# HGSOC, nPC = 10
# HNSCC, nPC = 20
# MEL, nPC = 19
# NSCLC early, nPC = 18
# NSCLC advanced, nPC = 30

### Cell clustering 
object <- FindNeighbors(object, dims = 1:nPC) 
res = c(0.4, 0.6, 0.8, 1, 1.2, 1.5)
object <- FindClusters(object, resolution = res)

### Run non-linear dimensional reduction (UMAP)
object <- RunUMAP(object, reduction = "pca", dims = 1:nPC)

### Plots
pdf("Resolutions.pdf")
for (i in 1:6){
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

# Selected resolution, with corresponding number of clusters, in each individual cancer type: 
# For BC_early, CRC, HCC, HNSCC and NSCLC_early, reference on corresponding publication (listed above)
# BC advanced: Resolution 0.4, 21 clusters
# CC: Resolution 0.5, 25 clusters
# GBM: Resolution 1, 34 clusters
# HGSOC: Resolution 0.2, 30 clusters
# MEL: Resolution 0.4, 31 clusters
# NSCLC advanced: Resolution 0.5, 21 clusters

### Find differentially expressed genes among clusters according to chosen resolution (e.g. 1) 
Idents(object) <- "RNA_snn_res.1"
DefaultAssay(object) <- "RNA"
object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=T)
write.table(object.markers, "object_markers.txt",col.names=NA, sep="\t")

###Plot marker genes to identify major cell types 
# Marker genes for major cell types shared among cancer types: 
# T-cell: "CD3D", "CD3E", "CD4", "CD8A"
# B-cell: "CD79A"
# Macrophages/monocytes: "CD68", "LYZ", "AIF1", "ITGAM"
# pDC: "LILRA4", "CXCR3"
# cDC: "CLEC9A", "XCR1", "CD1C", "CCR7", "CCL17", "CCL19"
# Mast cell: "MS4A2", "TPSAB1", "CPA3"
# Endothelial cell (EC): "CLDN5", "PECAM1", "VWF"
# Fibroblast: "COL1A1", "BGN", "DCN"
# Epithelial cell: "EPCAM", "KRT7", "KRT18"

# Cancer-specific markers have been used to identify cancer-specific cancer/epithelial cells and other cell types:
# enteric glia, erythroblasts, erythrocytes, muscle cells, oligodendrocytes, sinusoidal 

DefaultAssay(object) <- "RNA"
Idents(object) <- "RNA_snn_res.1"
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
Idents(object) <- "RNA_snn_res.1"
object <- RenameIdents(object,
                   `0` = "xx", `1` = "xx", `2` = "xx",`3` = "xx", `4` = "xx", `5` = "xx",
                   `6` = "xx", `7` = "xx", `8` = "xx",  `9` = "xx", `10` = "xx")
object$RNA_snn_res.1 <- Idents(object)
object$Majorcelltype_annotation <- object$RNA_snn_res.1

### Annotation visualization 
pdf("Annotation.pdf")
Idents(object) <- "Majorcelltype_annotation"
UMAPPlot(object)
dev.off()

### save Seurat Object with annotation
saveRDS(object, file = "/path/HGSOC.rds")
