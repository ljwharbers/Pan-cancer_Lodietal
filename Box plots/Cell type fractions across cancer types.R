### Load libraries 
library(Seurat)

#### Box plots displaying the fractions of major cell types detected across cancer types ####

###### Accessing high-quality single-cell data
# To generate a Seurat object containing all shared pan-cancer high-quality cells, users can access the read count data for each individual cancer type on our lab’s website (https://lambrechtslab.sites.vib.be/en/dataaccess).
# After performing annotation analysis (refer to "Individual Cancer Types Analysis" and the cell subtype annotations in the "Cell Subtype Annotation" folder), users can explore the metadata of our pan-cancer atlas, available as "Lodietall_metadata.csv" within the "Master_files" folder. 
# This metadata file includes both general sample information and detailed cell subtype annotations for further reference.
# Once this object has been generated, proceed with the following analysis. 

object <- readRDS("/path/object.RDS")

df <- object@meta.data

# how many samples? 
length(unique(dff$SampleID))
# 178

#### Keep only samples with more than 500 cells (see Methods)
dff = dff %>% group_by(SampleID) %>% mutate(n = n()) %>% filter(n >= 500)

# how many samples after selection? 
length(unique(dff$SampleID))
# 160

### Set the levels
# Compared to Majorcelltype_annotation, we further distinguish T-cells from NK-cells and separate Macrophages/Monocytes (Macro/Mono) from dendritic cells (DC).
# This refinement is achieved through sub typing analysis (see the files: "Tcell analysis", "Macrophage/Monocyte analysis" and "DC analysis").
# We refer to this intermediate level of annotation as "Intermediatecelltype_annotation"
my_levels <- c("T cell", "NK",  "B cell", "Macro/Mono", "DC", "Mast cell", "Cancer/epit cell", "Fibroblast","EC")
dff$Intermediatecelltype_annotation <- factor(x = dff$Intermediatecelltype_annotation, levels = my_levels)

p <- dff %>% 
  group_by(SampleID, TumorType, Intermediatecelltype_annotation) %>% 
  CalcFreq() %>%
  mutate(TumorType = factor(TumorType, levels = c("BC_early", "BC_adv", "CC", "CRC", 
                                                    "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv"))) %>%
  StatBoxPlot(ncol = 3, facet_scale = "free_y") + 
  labs(x = NULL, y=NULL)

celltype_column = "Intermediatecelltype_annotation"
data = p$data
#data[[celltype_column]] = factor(data[[celltype_column]])
groups = levels(data$TumorType)
celltypes = levels(data[[celltype_column]])
stat_test = data.frame(matrix(ncol = 12, nrow = 0))
colnames(stat_test) = c(celltype_column, ".y.", "group1", "group2", "n1", "n2", "statistic", 
                        "p", "y.position", "groups", "xmin", "xmax")

for (i in 1:length(groups)) {
  group1 = as.character(groups[i])
  for (j in i:length(groups)) {
    group2 = as.character(groups[j])
    
    if (group1 == group2)
      next
    
    # Calculate the test score
    for (ct in celltypes) {
      x = pull(data[data$TumorType == group1 & data[[celltype_column]] == ct, ], "freq")
      y = pull(data[data$TumorType == group2 & data[[celltype_column]] == ct, ], "freq")
      w = wilcox.test(x, y)
      
      row_ = data.frame(
        "celltype" = ct, 
        ".y." = "freq", 
        "group1" = group1, 
        "group2" = group2, 
        "n1" = length(x), 
        "n2" = length(y), 
        "statistic" = as.numeric(w$statistic), 
        "p" = w$p.value, 
        #"p.adj" = NA,
        #"p.adj.signif" = NA,
        "y.position" = 0, 
        "groups" = NA,
        "xmin" = min(i, j),
        "xmax" = max(i, j))
      #"p.signif" = NA,
      #"n" = NA)
      colnames(row_)[[1]] = celltype_column
      stat_test = rbind(stat_test, row_)
    }
  }
}

stat_test$Intermediatecelltype_annotation = factor(stat_test$Intermediatecelltype_annotation, levels = my_levels)

# Add significance
stat_test = stat_test %>% 
  filter(!is.na(stat_test$p)) %>%
  add_significance(p.col = "p", output.col = "p.signif") %>%
  # adjust_pvalue(p.col = "p", output.col = "p.adj", method="bonferroni") %>%
  adjust_pvalue(p.col = "p", output.col = "p.adj", method="holm") %>%
  # adjust_pvalue(p.col = "p", output.col = "p.adj", method= "none") %>%
  add_significance(p.col = "p.adj", output.col = "p.adj.signif") %>% 
  filter(p.adj.signif != "ns")

# Calculate statistic and show only top3 values
stat_test = stat_test %>% 
  group_by(!!sym(celltype_column)) %>%
  slice_min(p.adj, n=3, with_ties = F) %>% 
  ungroup()
max_vals = data %>% 
  group_by(!!sym(celltype_column)) %>% 
  slice_max(freq, n=1, with_ties = F) %>%
  dplyr::select(!!sym(celltype_column), freq)
stat_test = stat_test %>% 
  merge(max_vals, by=celltype_column) %>%
  group_by(!!sym(celltype_column)) %>% 
  mutate(n=seq_len(n())) %>%
  mutate(y.position = min(freq/10, 0.1)*n + freq) %>%
  ungroup()

options(repr.plot.width=15, repr.plot.height=8)
p = p + stat_pvalue_manual(data=stat_test, label = "p.adj.signif", tip.length = 0.01)

p$layers[[2]]$aes_params$size = 0.5
p$layers[[3]]$aes_params$size = 2

### Visualize the plot
p = p + scale_y_continuous(labels = percent_format(accuracy = 2), 
                           breaks = function(limits) {
                             x = 0.2
                             if (limits[2] < 0.4) {
                               x = 0.1
                             }
                             if (limits[2] < 0.2) {
                               x = 0.05
                             } 
                             if (limits[2] < 0.1) {
                               x = 0.02
                             }
                             return(seq(0, 1, x))
                           })
p = p + theme(strip.text.x = element_text(size=11))

options(repr.plot.width=20, repr.plot.height=12)
p

### Save the plot
ggsave("/path/Boxplot_intermediateltype_cellfraction_acrosscancertype.pdf", plot = p, width = 7, height = 6)




#### Box plots displaying the fractions of T-cell subtypes detected across cancer types ####
# As an example, we present the correlation of cell proportions across T-cell subtypes.
# The same approach has also been applied for the other subtypes analysis.

object <- readRDS("/path/object.RDS")

df <- object@meta.data

# how many samples? 
length(unique(dff$SampleID))
# 178

#### Keep only samples with more than 500 cells (see Methods)
dff = dff %>% group_by(SampleID) %>% mutate(n = n()) %>% filter(n >= 500)

# how many samples after selection? 
length(unique(dff$SampleID))
# 160

### Set the levels
# First, we add to the seurat object with all cells, the corresponding annotations at the cell sublusuter level (see the files: "Tcell analysis", "Bcell analysis",  "Macrophage/Monocyte analysis" and "DC analysis")
# We refer to this deep level of annotation as "Minorcelltype_annotation"
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
my_levels <- c(my_levels, levels(dff$Minorcelltype_annotation)[!(levels(dff$Minorcelltype_annotation) %in% my_levels)])
dff$Minorcelltype_annotation <- factor(x = dff$Minorcelltype_annotation, levels = my_levels)

### Selection of samples >= 20 T-cells (including NK) 
samples = dff %>% 
  group_by(SampleID, Majorcelltype_annotation) %>% 
  summarise(n=n()) %>% 
  filter(Majorcelltype_annotation == "T cell") %>% 
  filter(n >= 20) %>% 
  .$SampleID

# How many samples after selection?
dff = filter(dff, SampleID %in% samples)
# 153

p <- dff %>% 
  filter(Majorcelltype_annotation == "T cell") %>% droplevels() %>%
  group_by(SampleID, TumorType, Minorcelltype_annotation) %>% 
  CalcFreq() %>%
  mutate(TumorType = factor(TumorType, levels = c("BC_early", "BC_adv", "CC", "CRC", 
                                                    "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv"))) %>%
  StatBoxPlot(ncol = 6, facet_scale = "free_y") + 
  labs(x = NULL, y=NULL)

celltype_column = "Minorcelltype_annotation"
data = p$data
#data[[celltype_column]] = factor(data[[celltype_column]])
groups = levels(data$TumorType)
celltypes = levels(data[[celltype_column]])
stat_test = data.frame(matrix(ncol = 12, nrow = 0))
colnames(stat_test) = c(celltype_column, ".y.", "group1", "group2", "n1", "n2", "statistic", 
                        "p", "y.position", "groups", "xmin", "xmax")

for (i in 1:length(groups)) {
  group1 = as.character(groups[i])
  for (j in i:length(groups)) {
    group2 = as.character(groups[j])
    
    if (group1 == group2)
      next
    
    # Calculate the test score
    for (ct in celltypes) {
      x = pull(data[data$TumorType == group1 & data[[celltype_column]] == ct, ], "freq")
      y = pull(data[data$TumorType == group2 & data[[celltype_column]] == ct, ], "freq")
      w = wilcox.test(x, y)
      
      row_ = data.frame(
        "celltype" = ct, 
        ".y." = "freq", 
        "group1" = group1, 
        "group2" = group2, 
        "n1" = length(x), 
        "n2" = length(y), 
        "statistic" = as.numeric(w$statistic), 
        "p" = w$p.value, 
        #"p.adj" = NA,
        #"p.adj.signif" = NA,
        "y.position" = 0, 
        "groups" = NA,
        "xmin" = min(i, j),
        "xmax" = max(i, j))
      #"p.signif" = NA,
      #"n" = NA)
      colnames(row_)[[1]] = celltype_column
      stat_test = rbind(stat_test, row_)
    }
  }
}

stat_test$Minorcelltype_annotation = factor(stat_test$Minorcelltype_annotation, levels = my_levels)

# Add significance
stat_test = stat_test %>% 
  filter(!is.na(stat_test$p)) %>%
  add_significance(p.col = "p", output.col = "p.signif") %>%
  adjust_pvalue(p.col = "p", output.col = "p.adj", method="holm") %>%
  # adjust_pvalue(p.col = "p", output.col = "p.adj", method= "none") %>%
  add_significance(p.col = "p.adj", output.col = "p.adj.signif") %>% 
  filter(p.adj.signif != "ns")

# Calculate statistic and show only top3
stat_test = stat_test %>% 
  group_by(!!sym(celltype_column)) %>%
  slice_min(p.adj, n=3, with_ties = F) %>% 
  ungroup()
max_vals = data %>% 
  group_by(!!sym(celltype_column)) %>% 
  slice_max(freq, n=1, with_ties = F) %>%
  dplyr::select(!!sym(celltype_column), freq)
stat_test = stat_test %>% 
  merge(max_vals, by=celltype_column) %>%
  group_by(!!sym(celltype_column)) %>% 
  mutate(n=seq_len(n())) %>%
  mutate(y.position = min(freq/10, 0.1)*n + freq) %>%
  ungroup()

options(repr.plot.width=15, repr.plot.height=8)
p = p + stat_pvalue_manual(data=stat_test, label = "p.adj.signif", tip.length = 0.01)

p$layers[[2]]$aes_params$size = 0.5
p$layers[[3]]$aes_params$size = 2

### Visualize the plot
p = p + scale_y_continuous(labels = percent_format(accuracy = 2), 
                           breaks = function(limits) {
                             x = 0.2
                             if (limits[2] < 0.4) {
                               x = 0.1
                             }
                             if (limits[2] < 0.2) {
                               x = 0.05
                             } 
                             if (limits[2] < 0.1) {
                               x = 0.02
                             }
                             return(seq(0, 1, x))
                           })
p = p + theme(strip.text.x = element_text(size=11))

options(repr.plot.width=20, repr.plot.height=12)
p

### Save the plot
ggsave("/path/Boxplot_Tcellsubtypes_cellfraction_acrosscancertype.pdf", plot = p, width = 15, height = 6)
