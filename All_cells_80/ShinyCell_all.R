# Single cell RNA sequencing data in mouse melanoma tumors
library(Seurat)
library(ShinyCell)
library(dplyr)
library(ggsignif)
library(ggpubr)

setwd("~/OneDrive - KU Leuven/Projects/Shiny/pancancer_datasets/All_cells/")

# Load Takis data
seu = readRDS("../seurat_files/allcell_HQ_onlysharedcelltypes_pct40_sampled_seed111.rds")

seu@meta.data -> mt
colnames(mt)[c(3,15)]  <- c("technology_10x", "Cancer_type_early_adv")
mt$technology_10x <- as.factor(mt$technology_10x) 
mt$nCount_RNA <- NULL
mt$nFeature_RNA <- NULL
seu@meta.data <- mt


# Idents(seu) %>% unique() # 63 levels (Cell_subtype?)
# 
# prop.table(table(seu$Interm_cell_type))
# table(seu$Interm_cell_type)
# 
# # Subset per cell type - define first sample 
# celltypes <- unique(seu$Interm_cell_type)
# 
# list_objects <- list()
# Idents(seu) <- "Interm_cell_type"
# for (i in 1:length(celltypes)) {
#   
#   celltype <- celltypes[i]
#   print(celltype)
#   
#   object_celltype <- subset(seu, idents = celltype)
#   
#   n_cells <- Cells(object_celltype) %>% length()
#   
#   pct40 <- (n_cells * 40 / 100)
#   
#   set.seed(111)
#   sampled = sample(colnames(object_celltype), pct40)
#   
#   object_celltype_pct40 <- subset(object_celltype, cells = sampled)
#   
#   list_objects <- append(list_objects, object_celltype_pct40)
#   
# }
# 
# object_pct40_sampled <- merge(list_objects[[1]], list_objects[2:length(list_objects)])
# 
# prop.table(table(object_pct40_sampled$Interm_cell_type))
# table(object_pct40_sampled$Interm_cell_type)
# 
# seu_subset <- object_pct40_sampled
# rm(object_pct40_sampled)
## Need to generate PCA to be able to generate files
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu,npcs = 5)

dimensions <- 1:5
seu  <- FindNeighbors(seu , dims = dimensions)
seu  <- RunUMAP(seu,dims = dimensions)

DimPlot(seu,group.by = "Major_cell_type",label.box = T,label = T)
#Create config
scConf = createConfig(seu,maxLevels = 70)

#Define most relevant
scConf = modDefault(scConf, "Cell_subtype", "Cancer_type")



scConf = modColours(scConf, meta.to.mod = "Cancer_type",
                    new.colours= c("#c39beaff","#de5ea3ff","#bd5b60ff","#e49c6aff","#8fc887ff","#558cc6ff","#f8c352ff","#a0a1a2ff",
                                   "#40b2abff"))

scConf = modColours(scConf, meta.to.mod = "Cancer_type_early_adv",
                    new.colours= c("#9966c7ff","#d8b7f9ff","#de5ea3ff","#bd5b60ff","#e49c6aff","#8fc887ff","#558cc6ff","#f8c352ff","#a0a1a2ff",
                                   "#389b95ff","#84dcd2ff"))


scConf = modColours(scConf, meta.to.mod = "Cell_subtype",
                    new.colours= scales::hue_pal()(length(unique(seu$Cell_subtype))))


showLegend(scConf) 

#Remove lowly expressed genes
assay.data <- GetAssayData(object = seu,assay = "RNA")
above0 <- rowSums(assay.data > 0)
keep <- above0 > 5
keep <- keep[keep]
seu <- subset(seu, features = names(keep))

makeShinyFiles(seu, scConf, gene.mapping = TRUE,shiny.dir = "shinyApp/",shiny.prefix = "sc1",gex.assay = "RNA")

assay.data <- GetAssayData(object = seu,assay = "RNA")
s1 <- assay.data[,1:115000]
s2 <- assay.data[,115001:ncol(assay.data)]


saveRDS(assay.data,file = "shinyApp/assay_data.rds")
saveRDS(s1,file = "shinyApp/assay_data_sub1.rds")
saveRDS(s2,file = "shinyApp/assay_data_sub2.rds")


######### Deploy ##########
options(repos = BiocManager::repositories())
library(rsconnect)
options(rsconnect.max.bundle.size=5145728000)
deployApp("shinyApp/",appName = "PanCancer-scRNAseq")





