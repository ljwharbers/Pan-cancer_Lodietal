library(Seurat)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggdendro)
library(patchwork)
setwd("~/OneDrive - KU Leuven/Projects/Shiny/pancancer_datasets/T-cells/")
# Set randomizer seed for reproducibility
set.seed(1234)

figsize = function(width=8, height=8) {
  options(repr.plot.width=width, repr.plot.height=height)
}
filter_samples = TRUE

# Set folder where figures are saved
figure_folder = "figures_app"

pc2_original = readRDS("../PanCancer2_major_latest.rds")
pc2 = readRDS("../PanCancer2_major_latest.rds")
pc2

# Set the celltype column that should be used here
pc2$celltype = pc2$CellType_lev5

# For the purpose of this calculation NK cells are counted as part of the T-cell major group
pc2@meta.data[pc2$CellType_lev2.5 == "NK", "CellType_lev2.5"] = "T cell"

# Calculate the T-cell reactivity before subsetting
pc2 = AddModuleScore(pc2, features = list(c("ITGAE", "PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1", "TNFRSF9", "TNFRSF18")), name="Tcell_Reactivity")
pc2$Tcell_Reactivity = pc2$Tcell_Reactivity1
pc2$Tcell_Reactivity1 = NULL

#Filter samples with less than 500 cells
pc2$barcodes = colnames(pc2)
x = pc2@meta.data %>%
  filter(!(CellType_lev5 %in% c("Low quality", "Doublets", "Specific"))) %>%
  droplevels() %>%
  group_by(SampleID) %>% 
  mutate(n=n()) %>% 
  filter(n >= 500) %>%
  filter(BiopsySite %in% c("Tumor", "Metastasis"))

pc2 = subset(pc2, cells = x$barcodes)
pc2

colnames(pc2) -> keep_cells
saveRDS(file = "shinyApp/keep.rds" , object = keep_cells)

#Create list with all cell types
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

# ## Calculate again 
# # Calculate the T-cell reactivity before subsetting
pc2 = AddModuleScore(pc2, features = list(c("ITGAE", "PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1", "TNFRSF9", "TNFRSF18")), name="Tcell_Reactivity_after_filtering")
pc2$Tcell_Reactivity_filtering = pc2$Tcell_Reactivity_after_filtering1
pc2$Tcell_Reactivity_after_filtering1 = NULL


# Proliferative T-cells are also counted as T-cells when calculating the reactivity score; in other analyses they are sometimes excluded
df_react = pc2@meta.data %>% 
  filter(CellType_lev3 %in% c("T cell", "Proliferative T-cell")) %>% 
  group_by(SampleID) %>%
  summarise_at(.vars = vars(Tcell_Reactivity,Tcell_Reactivity_filtering), .funs = mean)

# Only keep samples with a reactivity score (so >0 Tcells)
df2_major = df2_major[rownames(df2_major) %in% df_react$SampleID, ]

# Merge the dataframes
df2 = merge(df2_major, df_react, by.x = 0, by.y = "SampleID") %>% 
  column_to_rownames("Row.names")

temp = pc2@meta.data
temp[temp$CellType_lev3 == "Proliferative T-cell", "CellType_lev3"] = "T cell"
tcell_count = temp %>% group_by(SampleID, CellType_lev3) %>% summarise(n=n()) %>% filter(CellType_lev3 == "T cell") %>% arrange(n)
too_few_tcells = tcell_count[tcell_count$n < 20, "SampleID"] %>% .$SampleID
df2 = df2[!rownames(df2) %in% too_few_tcells, ]

counts = x %>% 
  droplevels() %>%
  group_by(SampleID, CellType_lev2.5, .drop=F) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = "CellType_lev2.5", values_from = "n", values_fill = 0) %>% 
  column_to_rownames("SampleID")
colnames(counts) = paste0(colnames(counts), "_count")

df = merge(df2, counts, by=0) %>% column_to_rownames("Row.names")
df$SampleID = rownames(df)
df_copy = df

figsize(9, 8)
df_orig = df

# Keep only the samples with 20 or more T-cells (including the NK and proliferative T-cells)
df = if (filter_samples) filter(df_orig, `T cell_count` >= 20) else df_orig

cor_tcells = df[, c(tcell_types, "Tcell_Reactivity")] %>% 
  cor_test(vars = c(tcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
  filter(var2 == "Tcell_Reactivity") %>%
  filter(var1 != "Tcell_Reactivity")

# Keep only the samples with 10 or more B-cells
df = if (filter_samples) filter(df_orig, `B cell_count` >= 10) else df_orig

cor_bcells = df[, c(bcell_types, "Tcell_Reactivity")] %>% 
  cor_test(vars = c(bcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
  filter(var2 == "Tcell_Reactivity") %>%
  filter(var1 != "Tcell_Reactivity")

df = if (filter_samples) filter(df_orig, `Endothelial_count` >= 10) else df_orig

cor_endo = df[, c(endo_types[endo_types != "Specific"], "Tcell_Reactivity")] %>% 
  cor_test(vars = c(endo_types[endo_types != "Specific"], "Tcell_Reactivity"), method = "spearman") %>% 
  filter(var2 == "Tcell_Reactivity") %>%
  filter(var1 != "Tcell_Reactivity")

df = if (filter_samples) filter(df_orig, `Macrophages-Monocytes_count` >= 20) else df_orig

cor_myeloid = df[, c(macro_types, mono_types, "Neutrophils", "Tcell_Reactivity")] %>% 
  cor_test(vars = c(macro_types, mono_types, "Neutrophils", "Tcell_Reactivity"), method = "spearman") %>% 
  filter(var2 == "Tcell_Reactivity") %>%
  filter(var1 != "Tcell_Reactivity")

df = if (filter_samples) filter(df_orig, `DC_count` >= 10) else df_orig

cor_dc = df[, c(dc_types, "Tcell_Reactivity")] %>% 
  cor_test(vars = c(dc_types, "Tcell_Reactivity"), method = "spearman") %>% 
  filter(var2 == "Tcell_Reactivity") %>%
  filter(var1 != "Tcell_Reactivity")

ct_names = c(tcell_types, bcell_types, endo_types[endo_types != "Specific"], c(macro_types, mono_types, "Neutrophils"), dc_types)
major_types = c(
  rep("NK/T-cell", length(tcell_types)),
  rep("B-cell", length(bcell_types)),
  rep("Endothelial", length(endo_types[endo_types != "Specific"])),
  rep("Macrophages/Monocytes", length(c(macro_types, mono_types, "Neutrophils"))),
  rep("DC", length(dc_types))
)

# Create combined dataframe of all results
rvalues = c(cor_tcells$cor, cor_bcells$cor, cor_endo$cor, cor_myeloid$cor, cor_dc$cor)
pvalues = c(cor_tcells$p, cor_bcells$p, cor_endo$p, cor_myeloid$p, cor_dc$p)
df_2 = data.frame(celltype = ct_names, Score = rvalues, Pval = pvalues, name = "Spearman", CellType_lev2.5 = major_types)

figsize(30, 4)
df_2 = add_significance(df_2, p.col = "Pval", output.col = "sig")
df_2[df_2$sig == "ns", "sig"] = ""
df_2[df_2$sig == "****", "sig"] = "***"

pc_order = df_2 %>%
  arrange(desc(Score)) %>% 
  .$celltype
df_2$celltype = factor(df_2$celltype, levels=pc_order)
df_2$CellType_lev2.5 = factor(df_2$CellType_lev2.5, levels = c("NK/T-cell", "B-cell", "Macrophages/Monocytes", "DC", "Endothelial"))
df_2$name = "Pan-cancer"

# Draw the pancancer plot
p_pancancer = ggplot(df_2) + 
  geom_tile(aes(x = celltype, y = name, fill=Score), width=0.98, height=0.98) + 
  geom_text(aes(x = celltype, y = name, label = sig), size=5) + 
  facet_grid(~CellType_lev2.5, space = "free", scales = "free") + 
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-1, 1)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.text.y = element_blank(),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.x = element_text(size=14),
        strip.background = element_rect(fill="white", color="black"),
        axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5)) + 
  labs(x = NULL, y = NULL)

pdf(paste0(figure_folder, "/PC2_reactivity_filtered.pdf"), width=30, height=4)
p_pancancer + theme(axis.text.y = element_text(size=16))
dev.off()

p_pancancer + theme(axis.text.y = element_text(size=16))

res = NULL

ct_names =  c(tcell_types, bcell_types, endo_types[endo_types != "Specific"], c(macro_types, mono_types, "Neutrophils"), dc_types)
major_types = c(
  rep("NK/T-cell", length(tcell_types)),
  rep("B-cell", length(bcell_types)),
  rep("Endothelial", length(endo_types[endo_types != "Specific"])),
  rep("Macrophages/Monocytes", length(c(macro_types, mono_types, "Neutrophils"))),
  rep("DC", length(dc_types))
)

df_tumor = pc2@meta.data %>% 
  select(SampleID, TumorType3) %>% 
  distinct() %>% 
  filter(SampleID %in% rownames(df2))

# Iterate over the cancer types to calculate the correlations. Need to check if there are 
# any samples that pass the filtering and create empty dataframe if there's not.
for (cancer_type in unique(pc2$TumorType3)) {
  
  rvalues = c()
  pvalues = c()
  df_orig = df_copy
  
  # Keep only the samples of the selected tumor type
  samples = df_tumor %>% 
    filter(TumorType3 == cancer_type) %>%
    .$SampleID %>%
    as.vector()
  df_orig = df_orig[samples, ]
  
  # T-cells
  df = if (filter_samples) filter(df_orig, `T cell_count` >= 20) else df_orig
  if (nrow(df)) {
    cor_tcells = df[, c(tcell_types, "Tcell_Reactivity")] %>% 
      cor_test(vars = c(tcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
      filter(var2 == "Tcell_Reactivity") %>%
      filter(var1 != "Tcell_Reactivity")
  } else {
    cor_tcells = data.frame(cor = rep(NA, length(tcell_types)), 
                            p = rep(NA, length(tcell_types)))
  }
  
  # B-cells
  df = if (filter_samples) filter(df_orig, `B cell_count` >= 10) else df_orig
  if (nrow(df)) {
    cor_bcells = df[, c(bcell_types, "Tcell_Reactivity")] %>% 
      cor_test(vars = c(bcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
      filter(var2 == "Tcell_Reactivity") %>%
      filter(var1 != "Tcell_Reactivity")
  } else {
    cor_bcells = data.frame(cor = rep(NA, length(bcell_types)), 
                            p = rep(NA, length(bcell_types)))
  }
  
  # Endothelial
  df = if (filter_samples) filter(df_orig, `Endothelial_count` >= 10) else df_orig
  if (nrow(df)) {
    cor_endo = df[, c(endo_types[endo_types != "Specific"], "Tcell_Reactivity")] %>% 
      cor_test(vars = c(endo_types[endo_types != "Specific"], "Tcell_Reactivity"), method = "spearman") %>% 
      filter(var2 == "Tcell_Reactivity") %>%
      filter(var1 != "Tcell_Reactivity")
  } else {
    cor_endo = data.frame(cor = rep(NA, length(endo_types[endo_types != "Specific"])), 
                          p = rep(NA, length(endo_types[endo_types != "Specific"])))
  }
  
  # Myeloid
  df = if (filter_samples) filter(df_orig, `Macrophages-Monocytes_count` >= 20) else df_orig
  if (nrow(df)) {
    cor_myeloid = df[, c(macro_types, mono_types, "Neutrophils", "Tcell_Reactivity")] %>% 
      cor_test(vars = c(macro_types, mono_types, "Neutrophils", "Tcell_Reactivity"), method = "spearman") %>% 
      filter(var2 == "Tcell_Reactivity") %>%
      filter(var1 != "Tcell_Reactivity")
  } else {
    cor_myeloid = data.frame(cor = rep(NA, length(myeloid_types)), 
                             p = rep(NA, length(myeloid_types)))
  }
  
  # DC
  df = if (filter_samples) filter(df_orig, `DC_count` >= 10) else df_orig
  if (nrow(df)) {
    cor_dc = df[, c(dc_types, "Tcell_Reactivity")] %>% 
      cor_test(vars = c(dc_types, "Tcell_Reactivity"), method = "spearman") %>% 
      filter(var2 == "Tcell_Reactivity") %>%
      filter(var1 != "Tcell_Reactivity")
  } else {
    cor_dc = data.frame(cor = rep(NA, length(dc_types)), 
                        p = rep(NA, length(dc_types)))
  }
  
  # Create combined dataframe for the cancer type
  rvalues = c(cor_tcells$cor, cor_bcells$cor, cor_endo$cor, cor_myeloid$cor, cor_dc$cor)
  pvalues = c(cor_tcells$p, cor_bcells$p, cor_endo$p, cor_myeloid$p, cor_dc$p)
  df_2 = data.frame(celltype = ct_names, Score = rvalues, Pval = pvalues, name = "Spearman", CellType_lev2.5 = major_types)
  df_2 = add_significance(df_2, p.col = "Pval", output.col = "sig")
  df_2[df_2$sig == "ns", "sig"] = ""
  df_2[df_2$sig == "****", "sig"] = "***"
  temp = df_2
  temp$TumorType3 = cancer_type
  
  # Merge with large dataframe
  if(is.null(res))
    res = temp
  else 
    res = rbind(res, temp)
}

figsize(30, 9)

# Order by descending score of the pancancer data
res$celltype = factor(res$celltype, levels=pc_order)
res$CellType_lev2.5 = factor(res$CellType_lev2.5, levels = c("NK/T-cell", "B-cell", "Macrophages/Monocytes", "DC", "Endothelial"))
res$TumorType3 = factor(res$TumorType3, levels = rev(c("BC_early", "BC_adv", "CC", "CRC", "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv")))

# Draw the pancancer plot
p = ggplot(res) + 
  geom_tile(aes(x = celltype, y = TumorType3, fill = Score), width=0.98, height=0.98) + 
  geom_text(aes(x = celltype, y = TumorType3, label = sig), size=5) + 
  facet_grid(~CellType_lev2.5, space = "free", scales = "free") + 
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-1, 1)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.text.y = element_blank(),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.x = element_text(size=14),
        strip.background = element_rect(fill="white", color="black"),
        axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5)) + 
  labs(x = NULL, y = NULL)

p

p_pancancer = p_pancancer + theme(axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank()) + NoLegend()

figsize(30, 10)
p_temp = p + theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
p_temp = p_pancancer / p_temp + plot_layout(height=c(1, 7))
p_temp = p_temp & theme(axis.text.y = element_text(size=16))

pdf(paste0(figure_folder, "/PC2_reactivity_stratified_filtered.pdf"), width=30, height=10)
p_temp
dev.off()

df = p$data %>% 
  select(celltype, Score, TumorType3) %>%
  pivot_wider(names_from = "celltype", values_from = "Score", values_fill = 0) %>%
  column_to_rownames("TumorType3")

pc_df = p_pancancer$data
colnames(pc_df) = c("celltype", "Score", "Pval", "name", "CellType_lev2.5", "sig")
pc_df$TumorType3 = "Pan-cancer"

model <- hclust(dist(df))
dhc <- as.dendrogram(model)
# Rectangular lines

ddata <- dendro_data(dhc, type = "rectangle")

dp = ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) + 
  theme_void() +
  theme(panel.background = element_blank(), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 40,  # Bottom margin
                             l = 10))

p$data$TumorType3 = factor(p$data$TumorType3, levels = model$labels[model$order])

figsize(30, 10)
p_pc = p_pancancer + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(size=16)) + 
  NoLegend()

p_pc$data$name = "Pan-cancer"
p = p + theme(strip.text.x = element_blank(),
              strip.background = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
              axis.text.y = element_text(size=16))
fp = plot_spacer() + p_pc + dp + p + plot_layout(widths = c(1, 10), heights = c(1, 10))
fp = fp

pdf(paste0(figure_folder, "/PC2_reactivity_stratified_dendro_filtered.pdf"), width=30, height=10)
fp
dev.off()
