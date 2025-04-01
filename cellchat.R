######## Cellchat among cell subtypes belonging to TLS-like (analysis1) or type-1 immunity (analysis2) hub

### Load libraries 
library(Seurat)
library(CellChat)
library(patchwork)
library(openxlsx)
library(pheatmap)
options(stringsAsFactors = FALSE)

###### Accessing high-quality single-cell data
# To generate a Seurat object containing all shared pan-cancer high-quality cells, users can access the read count data for each individual cancer type on our lab’s website (https://lambrechtslab.sites.vib.be/en/dataaccess).
# After performing annotation analysis (refer to "Individual Cancer Types Analysis" and the cell subtype annotations in the "Cell Subtype Annotation" folder), users can explore the metadata of our pan-cancer atlas, available as "Lodietall_metadata.csv" within the "Master_files" folder. 
# This metadata file includes both general sample information and detailed cell subtype annotations for further reference.
# Once this object has been generated, proceed with the following analysis. 

#### Part I: Data input & processing and initialization of CellChat object ----

# Read object with all high-quality cells
object <- readRDS("/path/object.RDS")

# RDS-names to save at several steps
step_1 <- "Object_name_step1.rds"
step_2 <- "Object_name_step2.rds"
Final_step <- "Object_name_Final.rds"

# Analysis 1: Cellchat among cell subtypes of type-1 immunity gub
# my_levels <-  c("CD4+ TFH",	"Breg",	"Plasmablast",	"IFN Mac",	"IgG immature",	"IgG mature", "IgA immature", "IgA mature", "mQuiescDC", "GC B")
# Idents(Seurat_object) <- "Minorcelltype_annotation"
# Seurat_object <- subset(Seurat_object, idents = my_levels)
# Seurat_object$Minorcelltype_annotation <- factor(x = Seurat_object$Minorcelltype_annotation, levels = my_levels)

# Analysis 2: Cellchat among cell subtypes of TLS-like hub
my_levels <-  c("CD4+ TREG",	"mRegDC",	"Mono-like Mac",	"LAM2",	"Neutrophils",	"AXL_DC",	"Inflam Mac",	"CD4+ TH1",	"CD8+ TEX",	"Prolif T", "Lymphatic")
Idents(Seurat_object) <- "Minorcelltype_annotation"
Seurat_object <- subset(Seurat_object, idents = my_levels)
Seurat_object$Minorcelltype_annotation <- factor(x = Seurat_object$Minorcelltype_annotation, levels = my_levels)

Seurat_object <- subset(Seurat_object, subset = Identity %in% my_levels)
Seurat_object$CellChat_levels <- factor(x = Seurat_object$Identity, levels = my_levels)

Idents(Seurat_object) <- "CellChat_levels"

data.input <- GetAssayData(Seurat_object, assay = "RNA", layer = "data") 
labels <- Idents(Seurat_object)

# Create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels)) 

# Create a CellChat object using data matrix as input
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
print(levels(cellchat@idents))

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# Pre-processing the expression data for cell-cell communication analysis

# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)


##### Part II: Inference of cell-cell communication network ----
# Step 1
saveRDS(cellchat, paste0("/path/", step_1))

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, nboot = 10000)

# Step 2
saveRDS(cellchat, paste0("/path/", step_2))

# Filter out the cell-cell communication if there are only few number of cells (<10) in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Final step prior to plots
saveRDS(cellchat, paste0("/path/", Final_step))

# Check all the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways:
cellchat@netP[["pathways"]]



##### Part III: Plots ---- 

### Interaction weights/strenght and number of interactions (25%)
pdf("/path/network_interaction_plots_25pt.pdf", width = 10, height = 6)  
# Get group sizes  
groupSize <- as.numeric(table(cellchat@idents))  
# Set layout for two plots side by side  
par(mfrow = c(1,2), xpd=TRUE)  
options(repr.plot.width = 10, repr.plot.height = 6)  # Adjust width and height  
# Generate the plots  
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge= FALSE, title.name = "Interaction weights/strength (25%)", top=0.25)  
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge= FALSE, title.name = "Number of interactions (25%)", top=0.25)  
dev.off()   

### Chord plots for all singnificant pathways 
# Check the levels of cell identities (this will help to map the cell types correctly)
levels(cellchat@idents)

# Map cell types to broader groups: 
group.cellType <- c('CD4+ TFH'= "T",
                    'Breg'= "B",
                    'Plasmablast'= "B",
                    'IgG immature'= "B",
                    'IgG mature'= "B",
                    'IgA immature'= "B",
                    'IgA mature'= "B",
                    'GC B'= "B",
                    'mQuiescDC'= "DC",
                    'IFN Mac'= "Mac",
                    'PCV' = "EC")

# Ensure that all identities in cellchat@idents are covered by the group.cellType vector
group.cellType <- factor(group.cellType[levels(cellchat@idents)], 
                         levels = c("T", "B", "Mac", "DC", "EC"))

# Double-check the mapping of group.cellType
group.cellType

# Define the signaling pathways from: cellchat@netP[["pathways"]]
signaling_pathways <- c('MIF','MHC-II','APP','CypA','CD99','CXCL','GALECTIN','CD45','COLLAGEN','PECAM1',
                        'LAMININ','PECAM2','BAFF','CLEC','FN1','ICAM','MK','CCL','ANNEXIN','SELL','MHC-I',
                        'ApoE','NECTIN','SELE','IL16','CD40','SELPLG','Prostaglandin','ESAM','VISFATIN','LAIR1',
                        'SPP1','CDH5','JAM','CD34','SIRP','APRIL','THBS','CD86','BTLA','Cholesterol','ADGRE','TNF',
                        'PLAU','BAG','TENASCIN','ADGRG','SEMA4','CD23','CD6','PTPRM','GAP','PD-L1','NOTCH','CD46',
                        'PROS','GAS','MPZ','CADM')

# Open a single PDF file to save all plots
pdf("/path/chordplots_significantpathways_TLS-likehub.pdf", width = 7, height = 7)

# Loop through each signaling pathway and generate a chord plot
for (pathway in signaling_pathways) {
  par(mar = c(4, 4, 4, 4))  # Adjust margins to avoid clipping the title
  netVisual_chord_gene(cellchat, signaling = c(pathway), legend.pos.x = 5, legend.pos.y = 5)  # Adjust legend position
  title(main = pathway, cex.main = 1.5)  # Add pathway name as title
}

# Close the PDF device
dev.off()


