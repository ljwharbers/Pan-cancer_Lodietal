# We present the correlation among cell subtypes proportions and T-cell reactivity values, both pan-cancer and across cancer types.
# The same approach has also been applied for TEX and TREG analysis, as well as for the publicly available scRNAseq datasets collected and used to confirm our results.

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringi)
library(rstatix)
library(patchwork)
library(ggdendro)
library(textshape)

# Set randomizer seed for reproducibility (needed?)
set.seed(1234)

# Choose to filter samples
filter_samples = TRUE

# Specify directories
dir <- "/path/"
out_dir <- "/path/"

# For consistency, we created one file containing cell proportions ("Freq_SampleID.csv") and another containing T-cell reactivity scores ("TcellReactivity_SampleID.csv"), both calculated at the sample level.
# Specifically, in the first file, cell proportions are computed at the sample level:
# * For major cell types, proportions are calculated relative to the total number of cells.
# * At the cell subtype level, proportions are calculated relative to the sum of cells within the corresponding major cell type.
# In addition, we also have a file with the metadata of our cells ("metaData_ncells_SampleID.csv), in particular containing the total number of cells per sample.
# For more details, please refer to the files in the "Master files" folder.

# First, we refined cell type annotations through subclustering analysis. 
# Compared to Majorcelltype_annotation, we further distinguish T-cells from NK-cells and separate Macrophages/Monocytes (Macro/Mono) from DCs.
# This refinement is achieved through subclustering analysis (see the files: "Tcell analysis", "Macrophage/Monocyte analysis" and "DC analysis").
# We refer to this intermediate level of annotation as "Intermediatecelltype_annotation"

### Read data of the above-mentioned files with cell proportion
# Read file with T-cell reactivity scores per sample
signature_df <- read.csv(paste0(dir, "/master_files/Tcell_Reactivity_SampleID.csv"))

# Read file with cell proportions
freq_all_long <- read.csv(paste0(dir, "/master_files/Freq_SampleID.csv"))

ncells_T_NK <- read.csv(paste0(dir, "/master_files/ncells_T_NK_SampleID.csv"))
df_ncells_meta_T_NK <- dplyr::full_join(df_ncells_meta, ncells_T_NK)
df_ncells_meta_T_NK[is.na(df_ncells_meta_T_NK)] <- 0

# Read file with metadata 
df_ncells_meta <- read.csv(paste0(dir, "/master_files/metaData_ncells_SampleID.csv"))

# Filter frequencies file 
level = "Minorcelltype_annotation"
freq.long <- freq_all_long %>% filter(lev == level)

Minor_Intermediatecelltype_annotation <- read.csv(paste0(dir, "/master_files/Minor_Intermediatecelltype_annotation.csv"))
# Minor_2.0_Intermediatecelltype_annotation <- read.csv(paste0(dir, "/master_files/Minor_2.0_Intermediatecelltype_annotation.csv"))

# Fill in zero missing values
sample_freq <- freq.long %>%
  dplyr::select(-n, -lev) %>%
  tidyr::pivot_wider(names_from = "CellType", values_from = "freq", values_fill = NA) 

# Only keep samples with values (no NAs)
df_orig_all <- dplyr::left_join(signature_df, sample_freq)

# Add sample info
df_orig_all <- dplyr::left_join(df_orig_all, df_ncells_meta_T_NK)

### Filter frequencies based on annotation and samples (samples with >= 500 samples and only PrimaryTumor&Metastasis)
select = "500_T_M"
df_orig_all %>% dplyr::filter(ncells_noLQ_noD >= 500, BiopsySite %in% c("Tumor", "Metastasis")) %>% dim()
df_orig <- df_orig_all %>% dplyr::filter(ncells_noLQ_noD >= 500, BiopsySite %in% c("Tumor", "Metastasis"))

df_orig_CT <- df_orig %>% dplyr::filter(ncells_noLQ_noD_obj_MM >= 20)


### Define subtypes and celltypes 
TC_types = dplyr::filter(Minor_Intermediatecelltype_annotation, Intermediatecelltype_annotation == "TC") %>% dplyr::pull(Minorcelltype_annotation)
BC_types = dplyr::filter(Minor_Intermediatecelltype_annotation, Intermediatecelltype_annotation == "BC") %>% dplyr::pull(Minorcelltype_annotation)
MM_types = dplyr::filter(Minor_Intermediatecelltype_annotation, Intermediatecelltype_annotation == "MM") %>% dplyr::pull(Minorcelltype_annotation)
DC_types = dplyr::filter(Minor_Intermediatecelltype_annotation, Intermediatecelltype_annotation == "DC") %>% dplyr::pull(Minorcelltype_annotation)
EC_types = dplyr::filter(Minor_Intermediatecelltype_annotation, Intermediatecelltype_annotation == "EC") %>% dplyr::pull(Minorcelltype_annotation)
subtypes = c(TC_types, BC_types, MM_types, DC_types, EC_types)

major_types <- Minor_Intermediatecelltype_annotation %>% dplyr::rename(major_type = Intermediatecelltype_annotation, subtype = Minorcelltype_annotation)

major_types$major_type <- stri_replace_all_regex(major_types$major_type, pattern = c("TC", "BC", "MM", "DC", "EC"),
                                                 replacement = c("NK/T cells", "B cells",  "Monocytes and macrophages", "Dendritic cells", "Endothelial cells"), vectorize = FALSE)

### Correlation for all primary tumor and metastasis samples with >=500 cells per sample and >= 10 (B cells, DCs, ECs) or 20 cells (T cells and Monocytes and macrophages)
df = df_orig_CT
cor_subtypes = df[, c(subtypes, "avg_signature")] %>% 
  cor_test(vars = c(subtypes, "avg_signature"), method = "spearman") %>% 
  filter(var2 == "avg_signature") %>%
  filter(var1 != "avg_signature")

# Create combined dataframes of all results
cor_df = cor_subtypes %>% dplyr::rename(Score = cor, Pval = p, name = method, subtype = var1) %>% dplyr::left_join(major_types, by = "subtype") %>% dplyr::select(-var2)

# Go to "plot pancancer"


# Correlation for all primary tumor and metastasis samples with >=500 cells per sample and >= 10 (B cells, DCs, ECs) or 20 cells (T cells and Monocytes and macrophages)
select = "500_T_M_sub_filtered"

# T cell
# Keep only the samples with >= 20 or more T cells (excluding NKs)
df = if (filter_samples) filter(df_orig_CT, ncells_noLQ_noD_obj_TC >= 20) else df_orig_CT

cor_TC = df[, c(TC_types, "avg_signature")] %>% 
  cor_test(vars = c(TC_types, "avg_signature"), method = "spearman") %>% 
  filter(var2 == "avg_signature") %>%
  filter(var1 != "avg_signature")

# B cell
# Keep only the samples with >= 20 or more B cells
df = if (filter_samples) filter(df_orig_CT, ncells_noLQ_noD_obj_BC >= 10) else df_orig_CT

cor_BC = df[, c(BC_types, "avg_signature")] %>% 
  cor_test(vars = c(BC_types, "avg_signature"), method = "spearman") %>% 
  filter(var2 == "avg_signature") %>%
  filter(var1 != "avg_signature")

# Monocytes and macrophages
# Keep only the samples with 20 or more MMs 
df = if (filter_samples) filter(df_orig_CT, ncells_noLQ_noD_obj_MM >= 20) else df_orig_CT

cor_MM = df[, c(MM_types, "avg_signature")] %>% 
  cor_test(vars = c(MM_types, "avg_signature"), method = "spearman") %>% 
  filter(var2 == "avg_signature") %>%
  filter(var1 != "avg_signature")


# DC
# Keep only the samples with 10 or more DCs
df = if (filter_samples) filter(df_orig_CT, ncells_noLQ_noD_obj_DC >= 10) else df_orig_CT

cor_DC = df[, c(DC_types, "avg_signature")] %>% 
  cor_test(vars = c(DC_types, "avg_signature"), method = "spearman") %>% 
  filter(var2 == "avg_signature") %>%
  filter(var1 != "avg_signature")


# EC
# Keep only the samples with 10 or more ECs
df = if (filter_samples) filter(df_orig_CT, ncells_noLQ_noD_obj_EC >= 10) else df_orig_CT

cor_EC = df[, c(EC_types, "avg_signature")] %>% 
  cor_test(vars = c(EC_types, "avg_signature"), method = "spearman") %>% 
  filter(var2 == "avg_signature") %>%
  filter(var1 != "avg_signature")

# Create combined dataframe of all results
cor_subtypes <- bind_rows(list(cor_TC, cor_BC, cor_MM, cor_DC, cor_EC))
cor_df = cor_subtypes %>% dplyr::rename(Score = cor, Pval = p, name = method, subtype = var1) %>% dplyr::left_join(major_types, by = "subtype") %>% dplyr::select(-var2)

### Plot pan-cancer
cor_df = add_significance(cor_df, p.col = "Pval", output.col = "sig")
cor_df[cor_df$sig == "ns", "sig"] = ""
cor_df[cor_df$sig == "****", "sig"] = "***"

pc_order = cor_df %>%
  arrange(desc(Score)) %>% 
  .$subtype
cor_df$subtype = factor(cor_df$subtype, levels=pc_order)
cor_df$major_type = factor(cor_df$major_type, levels = c("NK/T cells", "B cells",  "Monocytes and macrophages", "Dendritic cells", "Endothelial cells"))
cor_df$name = "Pan-cancer"

### Draw the pan-cancer plot
p_pancancer = ggplot(cor_df) + 
  geom_tile(aes(x = subtype, y = name, fill=Score), width=0.98, height=0.98) + 
  geom_text(aes(x = subtype, y = name, label = sig), size=5) + 
  facet_grid(~major_type, space = "free", scales = "free") + 
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-1, 1)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.text.y = element_blank(),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.x = element_text(size=14),
        strip.background = element_rect(fill="white", color="black"),
        axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5)) + 
  labs(x = NULL, y = NULL)

### Save the plot
pdf(paste0(out_dir, "Correlation_heatmap_", signature, "_freq", level_index, "_", select, ".pdf"), width=20, height=3)
p_pancancer + theme(axis.text.y = element_text(size=16))  + ggtitle(signature) + labs(caption = select)
dev.off()

### Save correlation results 
write.csv(cor_df, paste0(out_dir, "Correlation_df_", signature, "_freq", level_index, "_", select, ".csv"), row.names = F)




### Stratification across cancer types 
# Correlation for all primary tumor and metastasis samples with >=500 cells per sample and >= 10 (B cells, DCs, ECs) or 20 cells (T cells and Monocytes and macrophages)
res = NULL

SampleID_TumorType <- df_orig_CT %>% dplyr::select(SampleID, TumorType) %>% unique

df_copy = df_orig_CT

for (cancer_type in unique(SampleID_TumorType$TumorType)) {
  
  # Keep only the samples of the selected tumor type
  samples = SampleID_TumorType %>% 
    dplyr::filter(TumorType == cancer_type) %>%
    .$SampleID %>%
    as.vector()
  
  df_orig_CT_cancerType = df_orig_CT %>% dplyr::filter(SampleID %in% samples)
  
  cor_subtypes = df_orig_CT_cancerType[, c(subtypes, "avg_signature")] %>% 
    cor_test(vars = c(subtypes, "avg_signature"), method = "spearman") %>% 
    filter(var2 == "avg_signature") %>%
    filter(var1 != "avg_signature")
  
  # Create combined dataframe for the cancer type
  cor_df = cor_subtypes %>% dplyr::rename(Score = cor, Pval = p, subtype = var1) %>% dplyr::left_join(major_types, by = "subtype") %>% dplyr::select(-var2)
  cor_df = add_significance(cor_df, p.col = "Pval", output.col = "sig")
  cor_df[cor_df$sig == "ns", "sig"] = ""
  cor_df[cor_df$sig == "****", "sig"] = "***"
  
  temp = cor_df
  temp$TumorType = cancer_type
  
  
  # Merge with large dataframe
  if(is.null(res))
    res = temp
  else 
    res = rbind(res, temp)
  
}

# Save correlation results 
write.csv(res, paste0(out_dir, "Correlation_df_", signature, "_freq", level_index, "_", select, "_stratifiedTumorType.csv"), row.names = F)

# Go to "plot pancancer and per tumortype"


# Correlation for all primary tumor and metastasis samples with >=500 cells per sample and >= 10 (B cells, DCs, ECs) or 20 cells (T cells and Monocytes and macrophages)
res = NULL

SampleID_TumorType <- df_ncells_meta %>% dplyr::select(SampleID, TumorType) %>% unique

df_copy = df_orig_CT

# Iterate over the cancer types to calculate the correlations
for (cancer_type in unique(SampleID_TumorType$TumorType)) {
  
  # Keep only the samples of the selected tumor type
  samples = SampleID_TumorType %>% 
    dplyr::filter(TumorType == cancer_type) %>%
    .$SampleID %>%
    as.vector()
  
  df_orig_CT_cancerType = df_orig_CT %>% dplyr::filter(SampleID %in% samples)
  
  
  # T cell
  df = if (filter_samples) filter(df_orig_CT_cancerType, ncells_noLQ_noD_obj_TC >= 20) else df_orig_CT_cancerType
  if (nrow(df)) {
    cor_TC = df[, c(TC_types, "avg_signature")] %>% 
      cor_test(vars = c(TC_types, "avg_signature"), method = "spearman") %>% 
      filter(var2 == "avg_signature") %>%
      filter(var1 != "avg_signature")
  } else {
    cor_TC = data.frame(var1 = TC_types,
                        var2 = "avg_signature",
                        cor = rep(NA, length(TC_types)),
                        statistic = rep(NA, length(TC_types)),
                        p = rep(NA, length(TC_types)),
                        method = "spearman")
  }
  
  
  # B cell
  df = if (filter_samples) filter(df_orig_CT_cancerType, ncells_noLQ_noD_obj_BC >= 10) else df_orig_CT_cancerType
  if (nrow(df)) {
    cor_BC = df[, c(BC_types, "avg_signature")] %>% 
      cor_test(vars = c(BC_types, "avg_signature"), method = "spearman") %>% 
      filter(var2 == "avg_signature") %>%
      filter(var1 != "avg_signature")
  } else {
    cor_BC = data.frame(var1 = BC_types,
                        var2 = "avg_signature",
                        cor = rep(NA, length(BC_types)),
                        statistic = rep(NA, length(BC_types)),
                        p = rep(NA, length(BC_types)),
                        method = "spearman")
  }
  
  
  # EC
  df = if (filter_samples) filter(df_orig_CT_cancerType, ncells_noLQ_noD_obj_EC >= 10) else df_orig_CT_cancerType
  if (nrow(df)) {
    cor_EC = df[, c(EC_types, "avg_signature")] %>% 
      cor_test(vars = c(EC_types, "avg_signature"), method = "spearman") %>% 
      filter(var2 == "avg_signature") %>%
      filter(var1 != "avg_signature")
  } else {
    cor_EC = data.frame(var1 = EC_types,
                        var2 = "avg_signature",
                        cor = rep(NA, length(EC_types)),
                        statistic = rep(NA, length(EC_types)),
                        p = rep(NA, length(EC_types)),
                        method = "spearman")
  }
  
  
  # Monocytes and macrophages
  df = if (filter_samples) filter(df_orig_CT_cancerType, ncells_noLQ_noD_obj_MM >= 20) else df_orig_CT_cancerType
  if (nrow(df)) {
    cor_MM = df[, c(MM_types, "avg_signature")] %>% 
      cor_test(vars = c(MM_types, "avg_signature"), method = "spearman") %>% 
      filter(var2 == "avg_signature") %>%
      filter(var1 != "avg_signature")
  } else {
    cor_MM = data.frame(var1 = MM_types,
                        var2 = "avg_signature",
                        cor = rep(NA, length(MM_types)),
                        statistic = rep(NA, length(MM_types)),
                        p = rep(NA, length(MM_types)),
                        method = "spearman")
  }
  
  
  # DC
  df = if (filter_samples) filter(df_orig_CT_cancerType, ncells_noLQ_noD_obj_DC >= 10) else df_orig_CT_cancerType
  if (nrow(df)) {
    cor_DC = df[, c(DC_types, "avg_signature")] %>% 
      cor_test(vars = c(DC_types, "avg_signature"), method = "spearman") %>% 
      filter(var2 == "avg_signature") %>%
      filter(var1 != "avg_signature")
  } else {
    cor_DC = data.frame(var1 = DC_types,
                         var2 = "avg_signature",
                         cor = rep(NA, length(DC_types)),
                         statistic = rep(NA, length(DC_types)),
                         p = rep(NA, length(DC_types)),
                         method = "spearman")
  }
  
  # Create combined dataframe for the cancer types
  cor_subtypes <- bind_rows(list(cor_TC, cor_BC, cor_MM, cor_DC, cor_EC))
  cor_df = cor_subtypes %>% dplyr::rename(Score = cor, Pval = p, subtype = var1) %>% dplyr::left_join(major_types, by = "subtype") %>% dplyr::select(-var2)
  
  cor_df = add_significance(cor_df, p.col = "Pval", output.col = "sig")
  cor_df[cor_df$sig == "ns", "sig"] = ""
  cor_df[cor_df$sig == "****", "sig"] = "***"
  
  temp = cor_df
  temp$TumorType = cancer_type
  
  # Merge with large dataframe
  if(is.null(res))
    res = temp
  else 
    res = rbind(res, temp)
}


# Save correlation results
write.csv(res, paste0(out_dir, "Correlation_df_", signature, "_freq", level_index, "_", select, "_stratifiedTumorType.csv"), row.names = F)



### plot pancancer and per tumortype

# Order by descending score of the pancancer data
res$subtype = factor(res$subtype, levels=pc_order)
res$major_type = factor(res$major_type, levels = c("NK/T cells", "B cells",  "Monocytes and macrophages", "Dendritic cells", "Endothelial cells"))
res$TumorType = factor(res$TumorType, levels = rev(c("BC_early", "BC_adv", "CC", "CRC", "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv")))

### Draw the pancancer plot
p = ggplot(res) + 
  geom_tile(aes(x = subtype, y = TumorType, fill = Score), width=0.98, height=0.98) + 
  geom_text(aes(x = subtype, y = TumorType, label = sig), size=5) + 
  facet_grid(~major_type, space = "free", scales = "free") + 
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
                                  axis.ticks.x = element_blank())

p_temp = p + theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
p_temp = p_pancancer / p_temp + plot_layout(height=c(1, 7))
p_temp = p_temp & theme(axis.text.y = element_text(size=16))

### Save the plot
pdf(paste0(out_dir, "Correlation_heatmap_", signature, "_freq", level_index, "_", select, "_stratifiedTumorType.pdf"), width=30, height=10)
p_temp
p_temp + ggtitle(signature) + labs(caption = select)
dev.off()

### Add dendrogram
df = p$data %>% 
  dplyr::select(subtype, Score, TumorType) %>%
  pivot_wider(names_from = "subtype", values_from = "Score", values_fill = 0) %>%
  tibble::column_to_rownames("TumorType")

pc_df = p_pancancer$data
pc_df$TumorType = "Pan-cancer"

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

p$data$TumorType = factor(p$data$TumorType, levels = model$labels[model$order])

p_pc = p_pancancer + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(size=16))

p_pc$data$name = "Pan-cancer"
p = p + theme(strip.text.x = element_blank(),
              strip.background = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
              axis.text.y = element_text(size=16))
fp = plot_spacer() + p_pc + dp + p + plot_layout(widths = c(1, 10), heights = c(1, 10))

### Save the plot
pdf(paste0(out_dir, "Correlation_heatmap_", signature, "_freq", level_index, "_", select, "_stratifiedTumorType_dendogram.pdf"), width=30, height=10)
fp
fp + ggtitle(signature) + labs(caption = select)
dev.off()
