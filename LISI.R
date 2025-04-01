# Info: https://github.com/immunogenomics/LISI

# Load libraries
library(lisi)
library(Seurat)
library(dplyr)
library(ggplot2)

###### Accessing high-quality single-cell data
# To generate Seurat object for each major cell type, users can access the read count data for each individual cancer type on our lab’s website (https://lambrechtslab.sites.vib.be/en/dataaccess).
# After performing annotation analysis (refer to  cell subtype annotations in the "Cell Subtype Annotation" folder), users can explore the metadata of our pan-cancer atlas, available as "Lodietall_metadata.csv" within the "Master_files" folder. 
# This metadata file includes both general sample information and detailed cell subtype annotations for further reference.
# Once this object has been generated, proceed with the following analysis. 

celltype <- "TCell"
# celltype <- "BCell"
# celltype <- "Macrophage/Monocyte"
# celltype <- "DC"
# celltype <- "EC"

object <- readRDS(paste0("/path/", celltype, ".rds"))

ann <- "Minorcelltype_annotation"


# Order levels in T-cells             
my_levels <- c("CD4+ TN", "CD4+ TEM", "CD4+ TFH", "CD4+ TH1", "CD4+ TREG", "CD4+ TH17", "CD8+ TN", "CD8+ TEM", "CD8+ TEX", "CD8+ TEMRA",
               "CD8+ TRM", "Prolif T", "CD8+ Tγδ", "CD8- Tγδ", "CD8+ MAIT", "MAIT", "NK inflam", "NK cyto")

# Order levels in B-cells
# my_levels <- c("Naive mature", "Memory IgM+", "Memory IgM-", "GC B", "Breg", "Plasmablast", "IgG mature","IgG immature", "IgA mature", "IgA immature")

# Order levels in Macrophage/Monocyte
# my_levels <- c("Classical Mono", "Non-classical Mono", "Mono-like Mac", "Inflam Mac", "IFN Mac", "Interm Mac", "Mono-like LAM",
# "LAM1", "Hypoxic Mac", "Perivasc Mac", "MT Mac", "LAM2", "Suppr Mac", "Prolif Mac",  "Neutrophils")

# Order levels in DCs
# my_levels <- c("cDC1", "cDC2", "AXL_DC", "DC3", "Lang-likeDC", "pDC", "mQuiescDC", "mRegDC")

# Order levels in ECs
# my_levels <- c("Arterial", "Capillary", "Stalk", "Tip", "Lymphatic", "IFN EC", "PCV", "Venous", "Prolif EC")

object$Minorcelltype_annotation <- factor(object$Minorcelltype_annotation, levels = my_levels)

# Plot
pdf(paste0("/path/UMAP_object_", celltype, ".pdf"), width = 5, height = 5)
DimPlot(object, group.by = ann, pt.size = 0.1, label.size = 5, raster = FALSE, repel = TRUE)
DimPlot(object, group.by = ann, label = TRUE, label.size = 5, pt.size = 0.1, raster = FALSE, repel = TRUE)  + theme_void() + ggtitle("") + theme(legend.position = "none")
DimPlot(object, group.by = ann, label = FALSE, label.size = 5, pt.size = 0.1, raster = FALSE, repel = TRUE) + theme_void() + ggtitle("") + theme(legend.position = "none") 
DimPlot(object, group.by = "Technology", label = FALSE, label.size = 5, pt.size = 0.1, raster = FALSE, repel = TRUE) + theme_void() + ggtitle("")
DimPlot(object, group.by = "Technology", label = FALSE, label.size = 5, pt.size = 0.1, raster = FALSE, repel = TRUE) + theme_void() + ggtitle("") + theme(legend.position = "none")
dev.off()

object$Technology <- as.factor(object$Technology)

# Meta data
object$barcode <- rownames(object@meta.data)
meta_data <- object[[]] %>% select(Technology)
rownames(meta_data) <- object$barcode

# Embeddings
X <- Embeddings(object = object, reduction = "umap")

lisi_res <- compute_lisi(X, meta_data, 'Technology')
colnames(lisi_res) <- "LISI"
lisi_res$barcode <- rownames(lisi_res)

# Join lisi
object@meta.data <- left_join(object@meta.data, lisi_res, by = 'barcode')

# Join umap embeddings
X <- data.frame(X)
X$barcode <- rownames(X)
object@meta.data <- left_join(object@meta.data, X, by = 'barcode')

rownames(object@meta.data) <- object$barcode


setwd("/path/LISI/")
pdf(paste0("LISI_", celltype, "_plot.pdf"), width = 5, height = 5)
ggplot(object@meta.data, aes(UMAP_1, UMAP_2, color = Technology)) +
  geom_point(shape = 21) #+facet_wrap(~key)

ggplot(object@meta.data, aes(UMAP_1, UMAP_2, color = LISI)) +
  geom_point(shape = 21) + ggtitle(celltype)+
  scale_fill_gradient(breaks=c(0,0.5,1, 1.5, 2),labels=c(0,0.5,1,1.5,2), limits=c(0,2))#+facet_wrap(~key)

ggplot(object@meta.data, aes(UMAP_1, UMAP_2, color = LISI)) +
  geom_point(shape = 21, size = 0.1) + ggtitle(celltype)+
  scale_fill_gradient(breaks=c(0,0.5,1, 1.5, 2),labels=c(0,0.5,1,1.5,2), limits=c(0,2)) + theme_void() + ggtitle("")

ggplot(object@meta.data, aes(UMAP_1, UMAP_2, color = LISI)) +
  geom_point(shape = 21, size = 0.1) + theme_void() + ggtitle("") + theme(legend.position = "none")

dev.off()


pdf(paste0("LISI_", celltype, "_density_plot.pdf"), width = 4, height = 2.5)
options(repr.plot.widt = 4, height = 2)
d <- density(object$LISI)
plot(d, xlim = c(0,2))
dev.off()

pdf(paste0("LISI_", celltype, "_boxplot.pdf"), width = 15, height = 7)
options(repr.plot.widt = 4, height = 2)
ggplot(object@meta.data, aes_string(x= ann, y = "LISI", fill = ann)) + geom_boxplot() + 
  NoLegend() + xlab("") + geom_jitter(alpha = 0.1, size = 0.5) + 
  #geom_point(alpha = 0.5, size = 1) +
  theme(text=element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
