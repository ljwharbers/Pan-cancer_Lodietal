# Single cell RNA sequencing data in mouse melanoma tumors
library(Seurat)
library(ShinyCell)
library(dplyr)
library(ggsignif)
library(ggpubr)

setwd("~/OneDrive - KU Leuven/Projects/Shiny/pancancer_datasets/BC/")

# Load Takis data
seu = readRDS("../seurat_files/bcell.rds")

seu@meta.data -> mt
colnames(mt)[c(3,15)]  <- c("technology_10x", "Cancer_type_early_adv")
mt$technology_10x <- as.factor(mt$technology_10x) 
mt$Major_cell_type <- NULL
mt$Interm_cell_type <- NULL
mt$Cell_subtype <- factor(mt$Cell_subtype, levels= levels(mt$Cell_subtype)[levels(mt$Cell_subtype) %in% unique(mt$Cell_subtype)])
mt$Tissue_type <- ifelse(mt$Tissue_type == T, "Tumor",'Normal')
seu@meta.data <- mt

#Create config
scConf = createConfig(seu)

#Define most relevant
scConf = modDefault(scConf, "Cell_subtype", "Cancer_type")


scConf = modColours(scConf, meta.to.mod = "Cancer_type",
                    new.colours= c("#c39beaff","#de5ea3ff","#bd5b60ff","#8fc887ff","#558cc6ff","#f8c352ff","#a0a1a2ff",
                                   "#40b2abff"))

scConf = modColours(scConf, meta.to.mod = "Cancer_type_early_adv",
                    new.colours= c("#9966c7ff","#d8b7f9ff","#de5ea3ff","#bd5b60ff","#8fc887ff","#558cc6ff","#f8c352ff","#a0a1a2ff",
                                   "#389b95ff","#84dcd2ff"))

scConf = modColours(scConf, meta.to.mod = "Cell_subtype",
                    new.colours= scales::hue_pal()(length(unique(seu$Cell_subtype))))


showLegend(scConf) 

# makeShinyApp(seu, scConf, gene.mapping = TRUE,#gex.assay = "SCT",
#              shiny.title = "PanCancer T-Cells")
# 
makeShinyFiles(seu, scConf, gene.mapping = TRUE,shiny.dir = "shinyApp/",shiny.prefix = "sc1",gex.assay = "RNA")

assay.data <- GetAssayData(object = seu,assay = "RNA")
assay.data <- Matrix(assay.data, sparse = TRUE)
saveRDS(assay.data,file = "shinyApp/assay_data.rds")

######### Deploy ##########
options(repos = BiocManager::repositories())
library(rsconnect)
options(rsconnect.max.bundle.size=5145728000)
deployApp("shinyApp/")
