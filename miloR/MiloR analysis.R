### Load libraries 
library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(Seurat)
library(ggpubr)
library(statmod)


### Set seed for reproducible results
set.seed(123)

### Define directory for output 
out_dir <- "/path/"


### Load object created through the analysis presented in "Tcell analysis", "Bcell analysis", "Macrophage/Monocyte analysis", "DC analysis", "EC analysis"
celltype <- "Tcell"
celltype <- "Bcell"
celltype <- "Macrophage/Monocyte"
celltype <- "DC"
celltype <- "EC"

object <- readRDS(paste0("/path/", celltype, ".rds"))


### Read metadata provided (provided in folder MiloR)
meta_Expansion <- read.csv("/path/metadata_MiloR.csv")

### Add expansion metadata to seurat object 
object$barcode <- rownames(object@meta.data)
object@meta.data <- dplyr::left_join(object@meta.data, meta_exp)
rownames(object@meta.data) <- object$barcode

### Make barplots about E (Expander) versus NE (Non-Expander)
x <- prop.table(table(object$Minorcelltype_annotation, object$Expansion), margin = 1)
x <- prop.table(table(object$Minorcelltype_annotation, object$Expansion), margin = 2)
df <- prop.table(table(object$Minorcelltype_annotation, object$Expansion), margin = 2) %>% as.data.frame()

df$Var1 <- factor(df$Var1, levels = my_levels)

pdf(paste0(out_dir, celltype, "_barplot_Expansion.pdf"), height = 5, width = 4)
ggplot(data=df, aes(x=Var2, y=Freq,fill=Var1)) +
  geom_bar(stat="identity") + 
  theme_classic() + theme(axis.text = element_text(color = "black"),
                          axis.title = element_text(color = "black"),
                          plot.title = element_text(color = "black"),
                          legend.text = element_text(color = "black"))
dev.off()


######## MiloR
### Create a Milo object ----
# In Bcell analysis, add:
# DefaultAssay(s) <- "integrated" 

### Convert to SingleCellExperiment
sce = as.SingleCellExperiment(object)

### Convert to Milo object
milo = Milo(sce)


### Construct KNN graph
# Set D to the number of PCA dimensions previously used in the Seurat analysis
# Tcell
K = 40 
D = 40 

# Bcell
K = 12
D = 12

# Macrophage/Monocyte
K = 17
D = 17

# DC
K = 11
D = 11

# EC
K = 13
D = 13

### buildGraph
# Bcell: 
# milo <- buildGraph(milo, k=K, d=D) 
milo <- buildGraph(milo, k=K, d=D, reduced.dim = "HARMONY")

### Constructing kNN graph
p=0.1
milo <- makeNhoods(milo, prop=p, k=K, d=D, refined=TRUE)

pdf(paste0(out_dir, celltype, "_plotNhoodSizeHist_", D, "D_", K, "K_", p, "p.pdf"))
plotNhoodSizeHist(milo)
dev.off()

milo$SampleID = factor(milo$SampleID)
milo <- countCells(milo, meta.data = object@meta.data, samples="SampleID") 

### Differential abundance testing 
design <- distinct(as.data.frame(milo@colData)[, c("SampleID", "Expansion", "TumorType")]) # add TumorType if used in design
colnames(design) = c("Sample", "Condition", "Tumortype")
rownames(design) = design$Sample

# Reorder rownames to match columns of nhoodCounts(milo)
design <- design[colnames(nhoodCounts(milo)), , drop=FALSE]
design %>% head

milo <- calcNhoodDistance(milo, d=D)

### Now we can do the test, explicitly defining our experimental design.
des <- "_TT1"

### testNhoods
da_results <- testNhoods(milo, design = ~ Tumortype + Condition, design.df = design) # add TT1 if used in design

### Using TMM normalisation
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

da_results %>%
  arrange(PValue) %>%
  head()

### Visualize neighbourhoods
# buildNhoodGraph
milo <- buildNhoodGraph(milo)

# plotNhoodGraphDA
pdf(paste0(out_dir, celltype, "_plotNhoodGraphDA_", D, "D_", K, "K_", p, "p", des, ".pdf"))
plotNhoodGraphDA(milo, milo_res=da_results, alpha=0.5) + scale_fill_gradient2(high = scales::muted("red"), mid = "white", low = scales::muted("blue"))  + labs(fill = "logFC")
dev.off()

# annotateNhoods
da_res <- annotateNhoods(milo, da_results, coldata_col = "Minorcelltype_annotation")
head(da_res)

# Plot fraction for final ann
pdf(paste0(out_dir, celltype, "_Minorcelltype_annotation_fraction_histogram_", D, "D_", K, "K_", p, "p", des, ".pdf"), height = 5, width = 5)
ggplot(da_res, aes(Minorcelltype_annotation_fraction)) + 
  geom_histogram(bins=50)
dev.off()

# Enforce extreme purity of clusters 
da_res$Minorcelltype_annotation_0.9 <- ifelse(da_res$Minorcelltype_annotation_fraction > 0.9, da_res$Minorcelltype_annotation, "Mixed")
da_res$Minorcelltype_annotation_0.8 <- ifelse(da_res$Minorcelltype_annotation_fraction > 0.8, da_res$Minorcelltype_annotation, "Mixed")
da_res$Minorcelltype_annotation_0.7 <- ifelse(da_res$Minorcelltype_annotation_fraction > 0.7, da_res$Minorcelltype_annotation, "Mixed")
head(da_res)

da_res$Minorcelltype_annotation <- factor(da_res$Minorcelltype_annotation, levels = rev(my_levels))
da_res$Minorcelltype_annotation_0.9 <- factor(da_res$Minorcelltype_annotation_0.9, levels = c("Mixed", rev(my_levels)))
da_res$Minorcelltype_annotation_0.8 <- factor(da_res$Minorcelltype_annotation_0.8, levels = c("Mixed", rev(my_levels)))
da_res$Minorcelltype_annotation_0.7 <- factor(da_res$Minorcelltype_annotation_0.7, levels = c("Mixed", rev(my_levels)))

# Plot DA beeswarm
pdf(paste0(out_dir, celltype, "_plotDAbeeswarm_", D, "D_", K, "K_", p, "p", des, ".pdf"))

plot = plotDAbeeswarm(da_res, group.by = "Minorcelltype_annotation", alpha = 0.1)
plot + scale_colour_gradient2(high = scales::muted("red"), mid = "white", low = scales::muted("blue")) + 
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"), 
        legend.text = element_text(color = "black")) # + xlab("")

plot = plotDAbeeswarm(da_res, group.by = "Minorcelltype_annotation_0.9", alpha = 0.1)
plot + scale_colour_gradient2(high = scales::muted("red"), mid = "white", low = scales::muted("blue")) + 
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"), 
        legend.text = element_text(color = "black")) # + xlab("")

plot = plotDAbeeswarm(da_res, group.by = "Minorcelltype_annotation_0.8", alpha = 0.1)
plot + scale_colour_gradient2(high = scales::muted("red"), mid = "white", low = scales::muted("blue")) + 
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"), 
        legend.text = element_text(color = "black")) # + xlab("")

plot = plotDAbeeswarm(da_res, group.by = "Minorcelltype_annotation_0.7", alpha = 0.1)
plot + scale_colour_gradient2(high = scales::muted("red"), mid = "white", low = scales::muted("blue")) + 
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"), 
        legend.text = element_text(color = "black")) # + xlab("")
dev.off()

### Save results
saveRDS(milo, paste0(out_dir, celltype, "_miloR_object_", D, "D_", K, "K_", p, "p", des, ".rds"))
write.table(da_res, paste0(out_dir, celltype, "_miloR_da_results_", D, "D_", K, "K_", p, "p", des, ".txt"), row.names = F)









