# Example on T-cells 

# Packages
```{r}
library(tidyr)
library(tidyverse)
library(rstatix)
library(Seurat)
library(RColorBrewer)

# Read object 
Tcells <- readRDS("/path/Tcells.rds")

# Calculate signature expression on T-cells 
naive <- readLines(con = "/path/naive.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(naive), ctrl = length(naive), name = 'naive') 

Activation_Effector_function <- readLines(con = "/path/Activation_Effector_function.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Activation_Effector_function), ctrl = length(Activation_Effector_function), name = 'Activation_Effector_function') 

Exhaustion <- readLines(con = "/path/Exhaustion.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Exhaustion), ctrl = length(Exhaustion), name = 'Exhaustion') 

TCR_signaling <- readLines(con = "/path/TCR_signaling.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(TCR_signaling), ctrl = length(TCR_signaling), name = 'TCR_signaling') 

Cytotoxicity <- readLines(con = "/Users/path/Cytotoxicity.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Cytotoxicity), ctrl = length(Cytotoxicity), name = 'Cytotoxicity') 

Cytokine_Cytokine_receptor <- readLines(con = "/path/Cytokine_Cytokine_receptor.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Cytokine_Cytokine_receptor), ctrl = length(Cytokine_Cytokine_receptor), name = 'Cytokine_Cytokine_receptor') 

Chemokine_Chemokine_receptor <- readLines(con = "/path/Chemokine_Chemokine_receptor.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Chemokine_Chemokine_receptor), ctrl = length(Chemokine_Chemokine_receptor), name = 'Chemokine_Chemokine_receptor') 

Senescence <- readLines(con = "/path/Senescence.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Senescence), ctrl = length(Senescence), name = 'Senescence') 

Anergy <- readLines(con = "/path/Anergy.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Anergy), ctrl = length(Anergy), name = 'Anergy') 

NFKB_signaling <- readLines(con = "/path/NFKB_signaling.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(NFKB_signaling), ctrl = length(NFKB_signaling), name = 'NFKB_signaling') 

Stress_response <- readLines(con = "/path/Stress_response.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Stress_response), ctrl = length(Stress_response), name = 'Stress_response') 

MAPK_signaling <- readLines(con = "/path/MAPK_signaling.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(MAPK_signaling), ctrl = length(MAPK_signaling), name = 'MAPK_signaling') 

Adhesion <- readLines(con = "/path/Adhesion.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Adhesion), ctrl = length(Adhesion), name = 'Adhesion') 

Costimulatory_molecules <- readLines(con = "/path/Costimulatory_molecules.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Costimulatory_molecules), ctrl = length(Costimulatory_molecules), name = 'Costimulatory_molecules') 

Treg_signature <- readLines(con = "/path/Treg_signature.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Treg_signature), ctrl = length(Treg_signature), name = 'Treg_signature') 

IFN_Response <- readLines(con = "/path/IFN_Response.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(IFN_Response), ctrl = length(IFN_Response), name = 'IFN_Response') 

Oxidative_phosphorylation <- readLines(con = "/path/Oxidative_phosphorylation.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Oxidative_phosphorylation), ctrl = length(Oxidative_phosphorylation), name = 'Oxidative_phosphorylation') 

Glycolysis <- readLines(con = "/path/Glycolysis.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Glycolysis), ctrl = length(Glycolysis), name = 'Glycolysis') 

Lipid_metabolism <- readLines(con = "/path/Lipid_metabolism.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Lipid_metabolism), ctrl = length(Lipid_metabolism), name = 'Lipid_metabolism') 

Fatty_acid_metabolism <- readLines(con = "/path/Fatty_acid_metabolism.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Fatty_acid_metabolism), ctrl = length(Fatty_acid_metabolism), name = 'Fatty_acid_metabolism') 

Pro_apoptosis <- readLines(con = "/path/Pro_apoptosis.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Pro_apoptosis), ctrl = length(Pro_apoptosis), name = 'Pro_apoptosis') 

Anti_apoptosis <- readLines(con = "/path/Anti_apoptosis.txt")
Tcells <- AddModuleScore(object = Tcells, features = list(Anti_apoptosis), ctrl = length(Anti_apoptosis), name = 'Anti_apoptosis') 



library(pheatmap)
colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
# add annotations 
features = c("naive1","Activation_Effector_function1", 
             "Exhaustion1", "TCR_signaling1", "Cytotoxicity1", "Cytokine_Cytokine_receptor1", "Chemokine_Chemokine_receptor1", "Senescence1",
             "Anergy1", "NFKB_signaling1", "Stress_response1", "MAPK_signaling1", "Adhesion1", "Costimulatory_molecules1", "Treg_signature1", "IFN_Response1", 
             "Oxidative_phosphorylation1", "Glycolysis1", "Lipid_metabolism1", "Fatty_acid_metabolism1", 
             "Pro_apoptosis1", "Anti_apoptosis1")

tf <- Tcells@meta.data[,c("annotation",features)] %>%
  group_by(annotation) %>% 
  summarise_all(~mean(.x, na.rm = TRUE)) %>%
  as.data.frame()

rownames(tf) <- tf$annotation
tf$annotation <- NULL

annotation = c("Differentiation", "Differentiation",
               "Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function","Function",
               "Metabolism","Metabolism","Metabolism","Metabolism",
               "Apoptosis","Apoptosis")

annotation = data.frame(Signatures = annotation)
rownames(annotation) = features

#specify colors 
Signatures = c("mediumseagreen", "sienna2", "mediumpurple3",  "darkgoldenrod2" )
names(Signatures) = c("Differentiation", "Function", "Metabolism",  "Apoptosis")
ann_colors = list(Signatures = Signatures)

#scaling for columns: normalize data within each column, which means to be able to compare cell types within the same column.
pheatmap(
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
  scale = "column", #scaling for columns
  angle_col = 90, # Rotate x-axis labels by 90 degrees
  labels_col = features  # Use modified features vector here
)
