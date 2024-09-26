# Box plots displaying the fractions of major cell types detected across cancer types. 
# Fractions were calculated for each sample with >500 cells (n=160) (see Methods). 

# Packages
```{r}
library(Seurat)
library(tidyr)
library(dplyr)

# Keep only high-quality and shared major cell types: read and filter Seurat object
Object <- readRDS("/path/all_cells.RDS)

Idents(Object) <- "major_annotation"
Object_selected <- subset(Object, idents = c("Doublets", "Low quality","Enteric glia", "Erythroblast", "Erythrocyte", "Muscle cell","Oligodendrocytes", "Sinusoidal"), invert = T)

# make a dataframe from the object
dff <- Object@meta.data

p <- dff %>% 
  group_by(SampleID, TumorType, major_annotation) %>% 
  CalcFreq() %>%
   mutate(TumorType3 = factor(TumorType, levels = c("BC_early", "BC_adv", "CC", "CRC", 
                                                    "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv"))) %>%
  StatBoxPlot(ncol = 3, facet_scale = "free_y") + 
  labs(x = NULL, y=NULL)

celltype_column = "major_annotation"
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
      x = pull(data[data$TumorType3 == group1 & data[[celltype_column]] == ct, ], "freq")
      y = pull(data[data$TumorType3 == group2 & data[[celltype_column]] == ct, ], "freq")
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

stat_test$major_annotation = factor(stat_test$major_annotation, levels = my_levels)

# Add significances
stat_test = stat_test %>% 
  filter(!is.na(stat_test$p)) %>%
  add_significance(p.col = "p", output.col = "p.signif") %>%
  adjust_pvalue(p.col = "p", output.col = "p.adj", method="holm") %>%
  add_significance(p.col = "p.adj", output.col = "p.adj.signif") %>% 
  filter(p.adj.signif != "ns")

# Calculate statistic and show only top 3. Update y.positions to get nice plots.
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

# Change some of the point sizes and theme
p$layers[[2]]$aes_params$size = 0.5
p$layers[[3]]$aes_params$size = 2
# Show plot
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
