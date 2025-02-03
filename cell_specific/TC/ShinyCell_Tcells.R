# Single cell RNA sequencing data in mouse melanoma tumors
library(Seurat)
library(ShinyCell)
library(dplyr)
library(ggsignif)
library(ggpubr)

setwd("~/OneDrive - KU Leuven/Projects/Shiny/pancancer_datasets/TC/")

# Load Takis data
seu = readRDS("../seurat_files/tcell.rds")
seu@reductions[["harmony"]] <- NULL

seu@meta.data -> mt
colnames(mt)[c(3,15)]  <- c("technology_10x", "Cancer_type_early_adv")
mt$technology_10x <- as.factor(mt$technology_10x) 
mt$Major_cell_type <- NULL
mt$Interm_cell_type <- as.character(mt$Interm_cell_type)  
mt$Cell_subtype <- factor(mt$Cell_subtype, levels= levels(mt$Cell_subtype)[levels(mt$Cell_subtype) %in% unique(mt$Cell_subtype)])
mt$Tissue_type <- ifelse(mt$Tissue_type == T, "Tumor",'Normal')
seu@meta.data <- mt

#Create config
scConf = createConfig(seu)

#Define most relevant
scConf = modDefault(scConf, "Cell_subtype", "Cancer_type")


scConf = modColours(scConf, meta.to.mod = "Cancer_type",
                    new.colours= c("#c39beaff","#de5ea3ff","#bd5b60ff","#e49c6aff","#8fc887ff","#558cc6ff","#f8c352ff","#a0a1a2ff",
                                   "#40b2abff"))

scConf = modColours(scConf, meta.to.mod = "Cancer_type_early_adv",
                    new.colours= c("#9966c7ff","#d8b7f9ff","#de5ea3ff","#bd5b60ff","#e49c6aff","#8fc887ff","#558cc6ff","#f8c352ff","#a0a1a2ff",
                                   "#389b95ff","#84dcd2ff"))

scConf = modColours(scConf, meta.to.mod = "Cell_subtype",
                    new.colours= c("#F8766D", "#E88526" ,"#D39200" ,"#B79F00", "#93AA00",
                                   "#5EB300" ,"#00BA38", "#00BF74", "#00C19F", "#00BFC4",
                                   "#00B9E3" ,"#00ADFA", "#619CFF", "#AE87FF", "#DB72FB",
                                   "#F564E3", "#FF61C3", "#FF699C"))
showLegend(scConf) 

# makeShinyApp(seu, scConf, gene.mapping = TRUE,#gex.assay = "SCT",
#              shiny.title = "PanCancer T-Cells")
# 
makeShinyFiles(seu, scConf, gene.mapping = TRUE,shiny.dir = "shinyApp/",shiny.prefix = "sc1")

  ################# Extra to make pre calculations for the plots ###########
library(tidyverse)
library("Hmisc")
library(corrplot)
library(patchwork)
set.seed(1234)

assay.data <- GetAssayData(object = seu)
saveRDS(assay.data,file = "shinyApp/assay_data.rds")

### Load seurat object with everything
pc2_orig = readRDS("../seurat_files/PanCancer2_major_latest.rds")
pc2 = pc2_orig
pc2_orig@meta.data$barcodes <- row.names(pc2_orig@meta.data)
saveRDS(file = "shinyApp/barcode.rds" , object = pc2_orig@meta.data[,c("barcodes","SampleID")])

## Use the global average to calculate controls
assay.data <- GetAssayData(object = pc2_orig)
data.avg <- Matrix::rowMeans(x = assay.data)
data.avg <- data.avg[order(data.avg)] ### Export this
saveRDS(data.avg,"shinyApp/data.avg.rds")


######################### Now generate files for correlation plots #########################
# Set the celltype column that should be used here
pc2$celltype = pc2$CellType_lev5
# For the purpose of this calculation NK cells are counted as part of the T-cell major group
pc2@meta.data[pc2$CellType_lev2.5 == "NK", "CellType_lev2.5"] = "T cell"

#Remove samples with less than 500 cells
pc2$barcodes = colnames(pc2)
x = pc2@meta.data %>%
  filter(!(CellType_lev5 %in% c("Low quality", "Doublets", "Specific"))) %>%
  droplevels() %>%
  group_by(SampleID) %>% 
  mutate(n=n()) %>% 
  filter(n >= 500) %>%
  filter(BiopsySite %in% c("Tumor", "Metastasis"))

pc2 = subset(pc2, cells = x$barcodes)

## Make a list with all main cell types and subypes
tcell_types = filter(pc2@meta.data, CellType_lev2 == "T cell") %>% 
  .$celltype %>% unique %>% sort
mono_types = filter(pc2@meta.data, CellType_lev4 == "Monocytes") %>% 
  .$celltype %>% unique %>% sort
macro_types = filter(pc2@meta.data, CellType_lev4 == "Macrophages") %>% 
  .$celltype %>% unique %>% sort
macro_types = c(macro_types, "Proliferating Macro")
myeloid_types = filter(pc2@meta.data, CellType_lev4 == "Monocytes") %>% 
  .$celltype %>% unique %>% sort
endo_types = filter(pc2@meta.data, CellType_lev2.5 == "Endothelial") %>% 
  .$celltype %>% unique %>% sort
bcell_types = filter(pc2@meta.data, CellType_lev2 == "B cell") %>% 
  .$celltype %>% unique %>% sort
dc_types = filter(pc2@meta.data, CellType_lev2.5 == "DC") %>% 
  .$celltype %>% unique %>% sort
celltypes = c(tcell_types, 
              bcell_types,
              myeloid_types,
              "Mast cell",
              dc_types)
other_types = unique(pc2$celltype)[!unique(pc2$celltype) %in% celltypes]
celltypes = c(celltypes, other_types)
celltypes = celltypes[!celltypes %in% "Specific"]

x = pc2@meta.data
celltypes_pc2 = unique(x$celltype)


all_types <- list("tcell_types" = tcell_types, "bcell_types" =  bcell_types,"myeloid_types" = myeloid_types,"Mast_cell" = c("Mast cell") ,
     "dc_types" = dc_types, "mono_types" = mono_types, "macro_types" = macro_types, "endo_types"  = endo_types )

# Normalize over major celltypes
df2_major = x %>% 
  group_by(SampleID, CellType_lev2.5, celltype) %>% 
  droplevels() %>%
  summarise(n=n()) %>%
  mutate(n=n/sum(n)) %>%
  ungroup() %>%
  select(-CellType_lev2.5) %>%
  as.data.frame() %>%
  filter(celltype %in% celltypes) %>%
  pivot_wider(names_from = "celltype", values_from = "n", values_fill = 0) %>%
  column_to_rownames("SampleID") %>%
  select(celltypes[!celltypes %in% c("Mast cell", "Cancer", "Fibroblast", "Epithelial")])

# Normalize over all cells for some columns that have no subtypes
df2_all = x %>%
  group_by(SampleID, celltype) %>% 
  droplevels() %>%
  summarise(n=n()) %>%
  mutate(n=n/sum(n)) %>%
  pivot_wider(names_from = "celltype", values_from = "n", values_fill = 0) %>%
  column_to_rownames("SampleID") %>%
  select(c("Mast cell", "Cancer", "Fibroblast", "Epithelial"))

# Merge the two dataframes
df2_major = merge(df2_major, df2_all, by=0) %>% column_to_rownames("Row.names")
df2_major = df2_major[rownames(df2_major) %in% sig1$SampleID, ]

temp = pc2@meta.data
temp[temp$CellType_lev3 == "Proliferative T-cell", "CellType_lev3"] = "T cell"
tcell_count = temp %>% group_by(SampleID, CellType_lev3) %>% summarise(n=n()) %>% filter(CellType_lev3 == "T cell") %>% arrange(n)
too_few_tcells = tcell_count[tcell_count$n < 20, "SampleID"] %>% .$SampleID

df2_major = df2_major[!rownames(df2_major) %in% too_few_tcells, ]

counts = x %>% 
  droplevels() %>%
  group_by(SampleID, CellType_lev2.5, .drop=F) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = "CellType_lev2.5", values_from = "n", values_fill = 0) %>% 
  column_to_rownames("SampleID")
colnames(counts) = paste0(colnames(counts), "_count")

df_tumor = pc2@meta.data %>% 
  select(SampleID, TumorType3) %>% 
  distinct() %>% 
  filter(SampleID %in% rownames(df2_major))

saveRDS(df2_major,file = "shinyApp/df2_major.rds")
saveRDS(df_tumor,file = "shinyApp/df_tumor.rds")
saveRDS(all_types,file = "shinyApp/alltypes.rds")
saveRDS(counts,file = "shinyApp/counts.rds")

library(msigdbr)
h_gene_sets = msigdbr(species = "Homo sapiens",category = "C2")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
saveRDS(msigdbr_list,file = "shinyApp/hallmarks.rds")
######### Deploy ##########
options(repos = BiocManager::repositories())
library(rsconnect)
options(rsconnect.max.bundle.size=5145728000)
deployApp("shinyApp/")
