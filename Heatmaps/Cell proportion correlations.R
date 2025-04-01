# Load libraries 
library(Hmisc)
library(ComplexHeatmap)
library(circlize)


#### Correlation of cell proportions across major cell types ####

###### Accessing high-quality single-cell data
# To generate a Seurat object containing all shared pan-cancer high-quality cells, users can access the read count data for each individual cancer type on our lab’s website (https://lambrechtslab.sites.vib.be/en/dataaccess).
# After performing annotation analysis (refer to "Individual Cancer Types Analysis" and the cell subtype annotations in the "Cell Subtype Annotation" folder), users can explore the metadata of our pan-cancer atlas, available as "Lodietall_metadata.csv" within the "Master_files" folder. 
# This metadata file includes both general sample information and detailed cell subtype annotations for further reference.
# Once this object has been generated, proceed with the following analysis. 

### Read object with all high-quality cells
object <- readRDS("/path/object.RDS")

df <- object@meta.data

# how many samples? 
length(unique(df$SampleID))
# 141 

### Keep only samples with more than 500 cells (see Methods)
df = df %>% group_by(SampleID) %>% mutate(n = n()) %>% filter(n >= 500)

# how many samples after selection? 
length(unique(df$SampleID))
# 125 

### Set the levels
# First, we refined cell type annotations through subclustering analysis. 
# Compared to Majorcelltype_annotation, we further distinguish T-cells from NK-cells and separate Macrophages/Monocytes (Macro/Mono) from DCs.
# This refinement is achieved through subclustering analysis (see the files: "Tcell analysis", "Macrophage/Monocyte analysis" and "DC analysis").
# We refer to this intermediate level of annotation as "Intermediatecelltype_annotation"

my_levels <- c("T cell", "NK",  "B cell", "Macro/Mono", "DC", "Mast cell", "Cancer/epit cell", "Fibroblast","EC")
df$Intermediatecelltype_annotation <- factor(x = df$Intermediatecelltype_annotation, levels = my_levels)

### Prepare to create the heatmap
Intermediatecelltype_annotation = df %>%
  droplevels() %>%
  group_by(SampleID, Intermediatecelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Intermediatecelltype_annotation", values_from = "freq", values_fill = 0)

# Reshape for use in correlation function
Intermediatecelltype_annotation[is.na(Intermediatecelltype_annotation)] = 0
Intermediatecelltype_annotation = column_to_rownames(Intermediatecelltype_annotation, "SampleID")

cols <- colorRampPalette(c("darkblue", "white", "darkred"))
cols <- cols(500)

# Create a dataframe of correlations for each pair of cell types. This will create a matrix of Spearman correlation coefficients between each pair of columns in the total_freq dataframe. 
corr_matrix <- rcorr(x = as.matrix(Intermediatecelltype_annotation), type = "spearman")
corr_r = corr_matrix$r
corr_p = corr_matrix$P

library(ComplexHeatmap)
library(circlize)

groups = list("T cell" =  c("CD4+ TN", "CD4+ TEM", "CD4+ TFH",
                            "CD4+ TH1", "CD4+ TREG", "CD4+ TH17",
                            "CD8+ TN", "CD8+ TEM", "CD8+ TEX", "CD8+ TEMRA",
                            "CD8+ TRM", "Prolif T",
                            "CD8+ Tγδ", "CD8- Tγδ", "CD8+ MAIT", "MAIT",
                            "NK inflam", "NK cyto"),
              "B cell" =  c("Naive mature", "Memory IgM+", "Memory IgM-", "GC B", "Breg", "Plasmablast", "IgG mature", 
                            "IgG immature", "IgA mature", "IgA immature"),
              "Macro/Mono" =  c("Classical Mono", "Non-classical Mono", "Mono-like Mac", "Inflam Mac", "IFN Mac", "Interm Mac", "Mono-like LAM",
                                "LAM1", "Hypoxic Mac", "Perivasc Mac", "MT Mac", "LAM2", "Suppr Mac", "Prolif Mac",  "Neutrophils"), 
              "DC" =  c("cDC1", "cDC2", "AXL_DC", "DC3", "Lang-likeDC", "pDC", "mQuiescDC", "mRegDC"),
              "EC" =  c("Arterial", "Capillary", "Stalk", "Tip", "Lymphatic", "IFN EC", "PCV", "Venous", "Prolif EC", " Specific_subtype"))

diag(corr_p) = 1

# Create the heatmap
heatmap <- Heatmap(
  corr_r,
  name = "Spearman Correlation",
  col = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")),
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = T,
  cluster_columns = T,
  rect_gp = gpar(col = "black", lwd = .1),
  column_title = "Correlations among major cell types",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(corr_p[i, j] < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    } else if(corr_p[i, j] < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    } else if(corr_p[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    }})

# Save the heatmap
ggsave("/path/Heatmap_correlation_cellproportions_Majorcelltypes.pdf", plot = heatmap) 


#### Correlation of cell proportions of T-cell subtypes versus all other subtypes (including T-cell subtypes) ####
# As an example, we present the correlation of cell proportions across T-cell subtypes.
# The same approach has also been applied for the other subtypes analysis.

# Read object with all cells (no doublets and no low quality cells)
object <- readRDS("/path/object.RDS")

df <- object@meta.data

# how many samples? 
length(unique(df$SampleID))
# 141 

# Keep only samples with more than 500 cells (see Methods)
df = df %>% group_by(SampleID) %>% mutate(n = n()) %>% filter(n >= 500)

# how many samples after selection? 
length(unique(df$SampleID))
# 125 

# Set the levels
my_levels <- c("Cancer/epithelial", "Fibroblast", "Mast cell",
               "CD4+ TN", "CD4+ TEM", "CD4+ TFH",
               "CD4+ TH1", "CD4+ TREG", "CD4+ TH17",
               "CD8+ TN", "CD8+ TEM", "CD8+ TEX", "CD8+ TEMRA",
               "CD8+ TRM", "Prolif T",
               "CD8+ Tγδ", "CD8- Tγδ", "CD8+ MAIT", "MAIT",
               "NK inflam", "NK cyto", 
               "Naive mature", "Memory IgM+", "Memory IgM-", "GC B", "Breg", "Plasmablast", "IgG mature", 
               "IgG immature", "IgA mature", "IgA immature",
               "Classical Mono", "Non-classical Mono", "Mono-like Mac", "Inflam Mac", "IFN Mac", "Interm Mac", "Mono-like LAM",
               "LAM1", "Hypoxic Mac", "Perivasc Mac", "MT Mac", "LAM2", "Suppr Mac", "Prolif Mac",  "Neutrophils",
               "cDC1", "cDC2", "AXL_DC", "DC3", "Lang-likeDC", "pDC", "mQuiescDC", "mRegDC",
               "Arterial", "Capillary", "Stalk", "Tip", "Lymphatic", "IFN EC", "PCV", "Venous", "Prolif EC")
df$Minorcelltype_annotation <- factor(x = df$Minorcelltype_annotation, levels = my_levels)

tcells_freq = df %>%
  filter(Majorcelltype_annotation == "T cell") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_MM = df %>%
  filter(Intermediatecelltype_annotation == "Macro/Mono") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_DC = df %>%
  filter(Intermediatecelltype_annotation == "DC") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_B = df %>%
  filter(Majorcelltype_annotation == "B cell") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_EC = df %>%
  filter(Majorcelltype_annotation == "EC") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

# I use all = TRUE to have a full join, not only the intersection of the 2 dataframes 
total_freq = merge(tcells_freq, extra_freq_MM, by="SampleID", all = TRUE)
total_freq = merge(total_freq, extra_freq_DC, by="SampleID", all = TRUE)
total_freq = merge(total_freq, extra_freq_B, by="SampleID", all = TRUE)
total_freq = merge(total_freq, extra_freq_EC, by="SampleID", all = TRUE)

# Reshape for use in correlation function. NA = 0 
total_freq[is.na(total_freq)] = 0
total_freq = column_to_rownames(total_freq, "SampleID")

cols <- colorRampPalette(c("darkblue", "white", "darkred"))
cols <- cols(500)

# create a dataframe of correlations for each pair of cell types 
corr_matrix <- rcorr(x = as.matrix(total_freq), type = "spearman")
corr_r = corr_matrix$r
corr_p = corr_matrix$P

groups = list("T cell" =  c("CD4+ TN", "CD4+ TEM", "CD4+ TFH",
                            "CD4+ TH1", "CD4+ TREG", "CD4+ TH17",
                            "CD8+ TN", "CD8+ TEM", "CD8+ TEX", "CD8+ TEMRA",
                            "CD8+ TRM", "Prolif T",
                            "CD8+ Tγδ", "CD8- Tγδ", "CD8+ MAIT", "MAIT",
                            "NK inflam", "NK cyto"),
              "B cell" =  c("Naive mature", "Memory IgM+", "Memory IgM-", "GC B", "Breg", "Plasmablast", "IgG mature", 
                            "IgG immature", "IgA mature", "IgA immature"),
              "Macro/Mono" =  c("Classical Mono", "Non-classical Mono", "Mono-like Mac", "Inflam Mac", "IFN Mac", "Interm Mac", "Mono-like LAM",
                                "LAM1", "Hypoxic Mac", "Perivasc Mac", "MT Mac", "LAM2", "Suppr Mac", "Prolif Mac",  "Neutrophils"), 
              "DC" =  c("cDC1", "cDC2", "AXL_DC", "DC3", "Lang-likeDC", "pDC", "mQuiescDC", "mRegDC"),
              "EC" =  c("Arterial", "Capillary", "Stalk", "Tip", "Lymphatic", "IFN EC", "PCV", "Venous", "Prolif EC"))

tcell_types = groups$`T cell`

diag(corr_p) = 1

subset_corr <- corr_r[, tcell_types]

row_split = c(rep("T cell", length(groups$`T cell`)),
              rep("Macro/Mono", length(groups$`Macro/Mono`)),
              rep("DC", length(groups$`DC`)),
              rep("B cell", length(groups$`B cell`)),
              rep("EC", length(groups$EC)))

# Create the heatmap
heatmap <- Heatmap(
  subset_corr,
  name = "Spearman Correlation",
  col = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_split = row_split,
  cluster_rows = T,
  cluster_columns = T,
  rect_gp = gpar(col = "black", lwd = .1),
  column_title = "Correlation of cell proportions of T-cell subtypes versus all other subtypes",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(corr_p[i, j] < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    } else if(corr_p[i, j] < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    } else if(corr_p[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    }})

# Save the heatmap
ggsave("/path/Heatmap_correlation_cellproportions_Tcellsversussubtypes.pdf", plot = heatmap) 





#### Correlation of cell proportions across all cell subtypes ####
# Read object with all cells (no doublets and no low quality cells)
object <- readRDS("/path/object.RDS")

df <- object@meta.data

# how many samples? 
length(unique(df$SampleID))
# 141 

# Keep only samples with more than 500 cells (see Methods)
df = df %>% group_by(SampleID) %>% mutate(n = n()) %>% filter(n >= 500)

# how many samples after selection? 
length(unique(df$SampleID))
# 125 

# Set the levels
my_levels <- c("Cancer/epithelial", "Fibroblast", "Mast cell",
               "CD4+ TN", "CD4+ TEM", "CD4+ TFH",
               "CD4+ TH1", "CD4+ TREG", "CD4+ TH17",
               "CD8+ TN", "CD8+ TEM", "CD8+ TEX", "CD8+ TEMRA",
               "CD8+ TRM", "Prolif T",
               "CD8+ Tγδ", "CD8- Tγδ", "CD8+ MAIT", "MAIT",
               "NK inflam", "NK cyto", 
               "Naive mature", "Memory IgM+", "Memory IgM-", "GC B", "Breg", "Plasmablast", "IgG mature", 
               "IgG immature", "IgA mature", "IgA immature",
               "Classical Mono", "Non-classical Mono", "Mono-like Mac", "Inflam Mac", "IFN Mac", "Interm Mac", "Mono-like LAM",
               "LAM1", "Hypoxic Mac", "Perivasc Mac", "MT Mac", "LAM2", "Suppr Mac", "Prolif Mac",  "Neutrophils",
               "cDC1", "cDC2", "AXL_DC", "DC3", "Lang-likeDC", "pDC", "mQuiescDC", "mRegDC",
               "Arterial", "Capillary", "Stalk", "Tip", "Lymphatic", "IFN EC", "PCV", "Venous", "Prolif EC")
df$Minorcelltype_annotation <- factor(x = df$Minorcelltype_annotation, levels = my_levels)

tcells_freq = df %>%
  filter(Majorcelltype_annotation == "T cell") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_MM = df %>%
  filter(Intermediatecelltype_annotation == "Macro/Mono") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_DC = df %>%
  filter(Intermediatecelltype_annotation == "DC") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_B = df %>%
  filter(Majorcelltype_annotation == "B cell") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

extra_freq_EC = df %>%
  filter(Majorcelltype_annotation == "EC") %>% # Keep only the cells you want to normalize for
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

# I use all = TRUE to have a full join, not only the intersection of the 2 dataframes 
total_freq = merge(tcells_freq, extra_freq_MM, by="SampleID", all = TRUE)
total_freq = merge(total_freq, extra_freq_DC, by="SampleID", all = TRUE)
total_freq = merge(total_freq, extra_freq_B, by="SampleID", all = TRUE)
total_freq = merge(total_freq, extra_freq_EC, by="SampleID", all = TRUE)

# Reshape for use in correlation function. NA = 0 
total_freq[is.na(total_freq)] = 0
total_freq = column_to_rownames(total_freq, "SampleID")


cols <- colorRampPalette(c("darkblue", "white", "darkred"))
cols <- cols(500)

corr_matrix <- rcorr(x = as.matrix(total_freq), type = "spearman")
corr_r = corr_matrix$r
corr_p = corr_matrix$P

diag(corr_p) = 1

# Create the heatmap
heatmap <- Heatmap(
  corr_r,
  name = "Spearman Correlation",
  col = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")),
  show_row_names = TRUE,
  show_column_names = TRUE,
  #row_split = row_split,
  cluster_rows = T,
  clustering_distance_rows = "euclidean", 
  clustering_method_rows = "complete",
  cluster_columns = T,
  clustering_distance_columns = "euclidean", 
  clustering_method_columns = "complete",
  rect_gp = gpar(col = "black", lwd = .1),
  column_title = "Correlation of cell proportions across all cell subtypes",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(corr_p[i, j] < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    } else if(corr_p[i, j] < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    } else if(corr_p[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 6),vjust = 0.75,hjust = 0.5)
    }})

# Save the heatmap
ggsave("/path/Heatmap_correlation_cellproportions_allsubtypes.pdf", plot = heatmap) 


