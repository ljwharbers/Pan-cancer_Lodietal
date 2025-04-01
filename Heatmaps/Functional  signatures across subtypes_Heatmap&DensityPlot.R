### Heatmap showing the (scaled) expression of curated functional gene signatures (columns) across subtypes (rows)

# As an example, we show the expression of functional gene signatures across T/NK-cell subtypes. The same approach has also been applied for the other subtypes analysis.

### Load libraries 
library(Seurat)
library(pheatmap)

### Read T/NK-cell object from individual cancer types (see file "Tcell analysis.R") 
Tcell <- readRDS("/path/Tcells.rds")

### Set the levels 
levels <- c("CD4+ TN", "CD4+ TEM", "CD4+ TFH",
             "CD4+ TH1", "CD4+ TREG", "CD4+ TH17",
             "CD8+ TN", "CD8+ TEM", "CD8+ TEX", "CD8+ TEMRA",
             "CD8+ TRM", "Prolif T",
             "CD8+ Tγδ", "CD8- Tγδ", "CD8+ MAIT", "MAIT",
             "NK inflam", "NK cyto")
Tcell$Minorcelltype_annotation <- factor(x = Tcell$Minorcelltype_annotation, levels = levels)

### First, create txt files from gene signatures from Supplemental Table 2 of Lodi et al.
### Then, read the files containing curated gene signatures and compute their average expression across cell subtypes
Naive <- readLines(con = "/path/Naive.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Naive), ctrl = length(Naive), name = 'Naive') 

Activation_Effector_function <- readLines(con = "/path/Activation_Effector_function.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Activation_Effector_function), ctrl = length(Activation_Effector_function), name = 'Activation_Effector_function') 

Exhaustion <- readLines(con = "/path/Exhaustion.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Exhaustion), ctrl = length(Exhaustion), name = 'Exhaustion') 

TCR_signaling <- readLines(con = "/path/TCR_signaling.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(TCR_signaling), ctrl = length(TCR_signaling), name = 'TCR_signaling') 

Cytotoxicity <- readLines(con = "/path/Cytotoxicity.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Cytotoxicity), ctrl = length(Cytotoxicity), name = 'Cytotoxicity') 

Cytokine_Cytokine_receptor <- readLines(con = "/path/Cytokine_Cytokine_receptor.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Cytokine_Cytokine_receptor), ctrl = length(Cytokine_Cytokine_receptor), name = 'Cytokine_Cytokine_receptor') 

Chemokine_Chemokine_receptor <- readLines(con = "/path/Chemokine_Chemokine_receptor.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Chemokine_Chemokine_receptor), ctrl = length(Chemokine_Chemokine_receptor), name = 'Chemokine_Chemokine_receptor') 

Senescence <- readLines(con = "/path/Senescence.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Senescence), ctrl = length(Senescence), name = 'Senescence') 

Anergy <- readLines(con = "/path/Anergy.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Anergy), ctrl = length(Anergy), name = 'Anergy') 

NFKB_signaling <- readLines(con = "/path/NFKB_signaling.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(NFKB_signaling), ctrl = length(NFKB_signaling), name = 'NFKB_signaling') 

Stress_response <- readLines(con = "/path/Stress_response.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Stress_response), ctrl = length(Stress_response), name = 'Stress_response') 

MAPK_signaling <- readLines(con = "/path/MAPK_signaling.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(MAPK_signaling), ctrl = length(MAPK_signaling), name = 'MAPK_signaling') 

Adhesion <- readLines(con = "/path/Adhesion.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Adhesion), ctrl = length(Adhesion), name = 'Adhesion') 

Costimulatory_molecules <- readLines(con = "/path/Costimulatory_molecules.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Costimulatory_molecules), ctrl = length(Costimulatory_molecules), name = 'Costimulatory_molecules') 

Treg_signature <- readLines(con = "/path/Treg_signature.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Treg_signature), ctrl = length(Treg_signature), name = 'Treg_signature') 

IFN_Response <- readLines(con = "/path/IFN_Response.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(IFN_Response), ctrl = length(IFN_Response), name = 'IFN_Response') 

Oxidative_phosphorylation <- readLines(con = "/path/Oxidative_phosphorylation.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Oxidative_phosphorylation), ctrl = length(Oxidative_phosphorylation), name = 'Oxidative_phosphorylation') 

Glycolysis <- readLines(con = "/path/Glycolysis.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Glycolysis), ctrl = length(Glycolysis), name = 'Glycolysis') 

Lipid_metabolism <- readLines(con = "/path/Lipid_metabolism.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Lipid_metabolism), ctrl = length(Lipid_metabolism), name = 'Lipid_metabolism') 

Fatty_acid_metabolism <- readLines(con = "/path/Fatty_acid_metabolism.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Fatty_acid_metabolism), ctrl = length(Fatty_acid_metabolism), name = 'Fatty_acid_metabolism') 

Pro_apoptosis <- readLines(con = "/path/Pro_apoptosis.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Pro_apoptosis), ctrl = length(Pro_apoptosis), name = 'Pro_apoptosis') 

Anti_apoptosis <- readLines(con = "/path/Anti_apoptosis.txt")
Tcell <- AddModuleScore(object = Tcell, features = list(Anti_apoptosis), ctrl = length(Anti_apoptosis), name = 'Anti_apoptosis') 

### Data scaling was performed to enable plot_density, as the function does not support negative values.
Tcell$Naive_score_scaled <-  round((Tcell$Naive1  + abs(min(Tcell$Naive1))) * 100)
Tcell$Activation_Effector_function_score_scaled <-  round((Tcell$Activation_Effector_function1  + abs(min(Tcell$Activation_Effector_function1))) * 100)
Tcell$Exhaustion_score_scaled <-  round((Tcell$Exhaustion1  + abs(min(Tcell$Exhaustion1))) * 100)
Tcell$TCR_signaling_score_scaled <-  round((Tcell$TCR_signaling1  + abs(min(Tcell$TCR_signaling1))) * 100)
Tcell$Cytotoxicity_score_scaled <-  round((Tcell$Cytotoxicity1  + abs(min(Tcell$Cytotoxicity1))) * 100)
Tcell$Cytokine_Cytokine_receptor_score_scaled <-  round((Tcell$Cytokine_Cytokine_receptor1  + abs(min(Tcell$Cytokine_Cytokine_receptor1))) * 100)
Tcell$Chemokine_Chemokine_receptor_score_scaled <-  round((Tcell$Chemokine_Chemokine_receptor1  + abs(min(Tcell$Chemokine_Chemokine_receptor1))) * 100)
Tcell$Senescence_score_scaled <-  round((Tcell$Senescence1  + abs(min(Tcell$Senescence1))) * 100)
Tcell$Anergy_score_scaled <-  round((Tcell$Anergy1  + abs(min(Tcell$Anergy1))) * 100)
Tcell$NFKB_signaling_score_scaled <-  round((Tcell$NFKB_signaling1  + abs(min(Tcell$NFKB_signaling1))) * 100)
Tcell$Stress_response_score_scaled <-  round((Tcell$Stress_response1  + abs(min(Tcell$Stress_response1))) * 100)
Tcell$MAPK_signaling_score_scaled <-  round((Tcell$MAPK_signaling1  + abs(min(Tcell$MAPK_signaling1))) * 100)
Tcell$Adhesion_score_scaled <-  round((Tcell$Adhesion1  + abs(min(Tcell$Adhesion1))) * 100)
Tcell$Costimulatory_molecules_score_scaled <-  round((Tcell$Costimulatory_molecules1  + abs(min(Tcell$Costimulatory_molecules1))) * 100)
Tcell$Treg_signature_score_scaled <-  round((Tcell$Treg_signature1  + abs(min(Tcell$Treg_signature1))) * 100)
Tcell$IFN_Response_score_scaled <-  round((Tcell$IFN_Response1  + abs(min(Tcell$IFN_Response1))) * 100)
Tcell$Oxidative_phosphorylation_score_scaled <-  round((Tcell$Oxidative_phosphorylation1  + abs(min(Tcell$Oxidative_phosphorylation1))) * 100)
Tcell$Glycolysis_score_scaled <-  round((Tcell$Glycolysis1  + abs(min(Tcell$Glycolysis1))) * 100)
Tcell$Lipid_metabolism_score_scaled <-  round((Tcell$Lipid_metabolism1  + abs(min(Tcell$Lipid_metabolism1))) * 100)
Tcell$Fatty_acid_metabolism_score_scaled <-  round((Tcell$Fatty_acid_metabolism1  + abs(min(Tcell$Fatty_acid_metabolism1))) * 100)
Tcell$Pro_apoptosis_score_scaled <-  round((Tcell$Pro_apoptosis1  + abs(min(Tcell$Pro_apoptosis1))) * 100)
Tcell$Anti_apoptosis_score_scaled <-  round((Tcell$Anti_apoptosis1  + abs(min(Tcell$Anti_apoptosis1))) * 100)

### Visualize with density plot  
pdf("Densityplot_Tcells_naive_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Naive_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Naive signature")
dev.off()      

pdf("Densityplot_Tcells_Activation_Effector_function_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Activation_Effector_function_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Activation/Effector function signature")
dev.off()      

pdf("Densityplot_Tcells_Exhaustion_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Exhaustion_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Exhaustion signature")
dev.off()      

pdf("Densityplot_Tcells_TCR_signaling_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "TCR_signaling_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("TCR signaling signature")
dev.off()      

pdf("Densityplot_Tcells_Cytotoxicity_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Cytotoxicity_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Cytotoxicity signature")
dev.off()      

pdf("Densityplot_Tcells_Cytokine_Cytokine_receptor_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Cytokine_Cytokine_receptor_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Cytokine/Cytokine receptor signature")
dev.off()      

pdf("Densityplot_Tcells_Chemokine_Chemokine_receptor_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Chemokine_Chemokine_receptor_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Chemokine/Chemokine receptor signature")
dev.off()      

pdf("Densityplot_Tcells_Senescence_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Senescence_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Senescence signature")
dev.off()      

pdf("Densityplot_Tcells_Anergy_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Anergy_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Anergy signature")
dev.off()      

pdf("Densityplot_Tcells_NFKB_signaling_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "NFKB_signaling_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("NFKB signaling signature")
dev.off()      

pdf("Densityplot_Tcells_Stress_response_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Stress_response_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Stress response signature")
dev.off()      

pdf("Densityplot_Tcells_MAPK_signaling_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "MAPK_signaling_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("MAPK signaling signature")
dev.off()      

pdf("Densityplot_Tcells_Adhesion_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Adhesion_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Adhesion signature")
dev.off()      

pdf("Densityplot_Tcells_Costimulatory_molecules_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Costimulatory_molecules_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Costimulatory molecules signature")
dev.off()      

pdf("Densityplot_Tcells_Treg_signature_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Treg_signature_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Treg signature")
dev.off()      

pdf("Densityplot_Tcells_IFN_Response_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "IFN_Response_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("IFN Response signature")
dev.off()      

pdf("Densityplot_Tcells_Oxidative_phosphorylation_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Oxidative_phosphorylation_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Oxidative phosphorylation signature")
dev.off()      

pdf("Densityplot_Tcells_Glycolysis_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Glycolysis_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Glycolysis signature")
dev.off()      

pdf("Densityplot_Tcells_Lipid_metabolism_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Lipid_metabolism_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Lipid metabolism signature")
dev.off()      

pdf("Densityplot_Tcells_Fatty_acid_metabolism_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Fatty_acid_metabolism_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Fatty acid metabolism signature")
dev.off()      

pdf("Densityplot_Tcells_Pro_apoptosis_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Pro_apoptosis_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Pro apoptosis signature")
dev.off()      

pdf("Densityplot_Tcells_Anti_apoptosis_magma.pdf", width = 10, height = 9)          
plot_density(Tcell, "Anti_apoptosis_score_norm", reduction = "umap", pal = "magma", method = "wkde") +ggtitle("Anti apoptosis signature")
dev.off()      

### Prepare to create the heatmap
# Set the color palette
colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

### Add annotations 
features = c("Naive1","Activation_Effector_function1", 
             "Exhaustion1", "TCR_signaling1", "Cytotoxicity1", "Cytokine_Cytokine_receptor1", "Chemokine_Chemokine_receptor1", "Senescence1",
             "Anergy1", "NFKB_signaling1", "Stress_response1", "MAPK_signaling1", "Adhesion1", "Costimulatory_molecules1", "Treg_signature1", "IFN_Response1", 
             "Oxidative_phosphorylation1", "Glycolysis1", "Lipid_metabolism1", "Fatty_acid_metabolism1", 
             "Pro_apoptosis1", "Anti_apoptosis1")

tf <- Tcell@meta.data[,c("Minorcelltype_annotation",features)] %>%
  group_by(Minorcelltype_annotation) %>% 
  summarise_all(~mean(.x, na.rm = TRUE)) %>%
  as.data.frame()

rownames(tf) <- tf$Minorcelltype_annotation
tf$Minorcelltype_annotation <- NULL

annotation = c("Differentiation", "Differentiation",
               "Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function",
               "Metabolism","Metabolism","Metabolism","Metabolism",
               "Apoptosis","Apoptosis")

annotation = data.frame(Signatures = annotation)
rownames(annotation) = features

### Specify colors per set of signatures 
Signatures = c("mediumseagreen", "sienna2", "mediumpurple3",  "darkgoldenrod2" )
names(Signatures) = c("Differentiation", "Function", "Metabolism",  "Apoptosis")
ann_colors = list(Signatures = Signatures)

### Create the heatmap
heatmap <- pheatmap(
  tf,
  color = colors,
  treeheight_col = 0,
  cluster_cols = FALSE,
  border_color = TRUE,
  cellwidth = 15,
  cellheight = 15,
  fontsize = 12,
  main = "T-cells functional signatures",
  annotation = annotation, 
  annotation_colors = list(Signatures = Signatures),
  scale = "column", # Scaling for columns
  angle_col = 90, # Rotate x-axis labels by 90 degrees
  labels_col = features  # Use modified features vector
)

### Save the heatmap
ggsave("/path/Heatmap_signatures_Tcells.pdf", plot = heatmap) 
