# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(forcats)
library(tibble)

# Specify directories
dir <- "/path/"
out_dir <- "/path/"

# For consistency, we created one file containing cell proportions ("Freq_SampleID.csv") and another containing Tcell expansion data in early BC and HNSCC ("meta_object_expansion_only_T_M_noNA.csv"), both calculated at the sample level.
# Specifically, in the first file, cell proportions are computed at the sample level:
# * For major cell types, proportions are calculated relative to the total number of cells.
# * At the cell subtype level, proportions are calculated relative to the sum of cells within the corresponding major cell type.
# In addition, we also have a file with the metadata of our cells ("metaData_ncells_SampleID.csv) with the total number of cells per sample. For more details, please refer to the files in the "Master_files" folder.

### Read data of the above-mentioned files 

# Read expansion data
meta_object_expansion <- read.csv(paste0(dir, "/Master_files/meta_object_expansion_only_T_M_noNA.csv")) 

# Read file with metadata 
df_ncells_meta <- read.csv(paste0(dir, "/Master_files/metaData_ncells_SampleID.csv"))

# Read file with cell proportions
freq_all_long <- read.csv(paste0(dir, "/Master_files/freq_all_SampleID_long.csv"))

# Filter cell frequencies file based on expansion
level = "Minorcelltype_annotation"
freq.long <- freq_all_long %>% filter(lev == level)
freq.long <- freq.long %>% dplyr::filter(SampleID %in% samples)

#### Select cell subtypes belonging to Type-1 immunity hub
freq.long <- freq.long %>% dplyr::filter(CellType %in%
                                           c("CD4+ TREG", "mRegDC", "Mono-like Mac", "LAM2", "Neutrophils", "AXL_DC", "Inflam Mac", "Lymphatic",
                                             "CD4+ TH1", "CD8+ TEX", "Prolif T"))
hub <- "Type-1 immunity hub"

### Select cell subtypes belonging to TLS-like hub
# freq.long <- freq.long %>% dplyr::filter(CellType %in%
#                                            c("CD4+ TFH", "GC B","Breg", "Plasmablast", "IFN Mac",
#                                              "IgG mature", "IgG immature", "IgA mature", "IgA immature", "mQuiescDC"))
# 
# 
# hub <- "TLS-like hub"

### Add missing values
freq.long0 <- freq.long %>% tidyr::complete(SampleID, CellType, lev)
freq.long0[is.na(freq.long0)] <- 0

### Combine cell subtypes of the corresponding hub
sample_freq_sum <- freq.long0 %>% dplyr::group_by(SampleID) %>% dplyr::summarise(n_sum = sum(n), freq_sum = sum(freq))

ncells <- df_ncells_meta %>% dplyr::select(SampleID, ncells, ncells_noLQ_noD)

sample_freq_sum <- dplyr::left_join(sample_freq_sum, ncells)
sample_freq_sum$combined_prop.clean <- sample_freq_sum$n_sum / sample_freq_sum$ncells_noLQ_noD

### Merge proportions with expansion data (filtered already for expansion only)
sample_freq_sum <- dplyr::full_join(sample_freq_sum, meta_object_expansion, by = "SampleID") 
included <- "all_T_M"

df <- sample_freq_sum
select <- paste0(hub, "_all_T_M")


### Select only samples with >= 500 samples and only PrimaryTumor&Metastasis 
select <- paste0(hub, "_500_T_M")
df <-  sample_freq_sum %>% filter(ncells >= 500) 
df <-  sample_freq_sum %>% filter(ncells_noLQ_noD >= 500) 

df$expansion <- factor(df$expansion, levels = c("NE", "E"))

### Create and save the plot
colors <- c("midnightblue", "red2")
y_axis <- "combined_prop.clean"

pdf(paste0(out_dir, "patient_boxplot_", select, "_expansion.pdf"), width = 4, height = 4)
 plot <- ggplot(df, aes_string(x= "expansion", y= y_axis, color = "expansion")) + 
    scale_y_continuous(expand = expansion(mult = c(.01, .1))) +
    theme_classic() +
    xlab("") + ggtitle("") +
    scale_color_manual(name = "expansion", values = colors) +
    theme(text=element_text(size=16, colour = "black")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.ticks = element_blank(), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    geom_boxplot(position=position_dodge(1), fill = NA) + 
    geom_jitter(shape=16, position = position_dodge(1))
  
  print(plot + stat_compare_means(label = "p.format", method = "wilcox.test", paired = F, show.legend = F, size = 5) + labs(caption = select))
  print(plot + stat_compare_means(label = "p.signif", method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 5, label.y = 0.95*max(df[,y])) + theme(text=element_text(size=16)) + labs(caption = select))
  print(plot + stat_compare_means(label = "p.signif", method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 5) + theme(text=element_text(size=16)) + labs(caption = select))

  dev.off()

#### Across cancer type
pdf(paste0(out_dir, "patient_boxplot_", select, "_TumorType_expansion.pdf"), width = 6, height = 4)

plot <- ggplot(df, aes_string(x= "expansion", y= y_axis, color = "expansion")) + 
    scale_y_continuous(expand = expansion(mult = c(.01, .1))) +
    theme_classic() +
    xlab("") + ggtitle("") +
    scale_color_manual(name = "expansion", values = colors) +
    theme(text=element_text(size=16, colour = "black")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.ticks = element_blank(), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    geom_boxplot(position=position_dodge(1), fill = NA) + 
    geom_jitter(shape=16, position = position_dodge(1)) +
    facet_grid(~TumorType, scales = "free", space = "free_x",)
  
  print(plot + stat_compare_means(label = "p.format", method = "wilcox.test", paired = F, show.legend = F, size = 5) + labs(caption = select))
  print(plot + stat_compare_means(label = "p.signif", method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 5, label.y = 0.95*max(df[,y])) + theme(text=element_text(size=16)) + labs(caption = select))
  print(plot + stat_compare_means(label = "p.signif", method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 5) + theme(text=element_text(size=16)) + labs(caption = select))

  dev.off()

  
#### Scatter plots  
cols <- c(
  HNSCC = "#F2B342",
  BC_early = "#E0CCEE"
)

pdf(paste0(out_dir, "correlation_nExp2_", select, ".pdf"), width = 5, height = 5)
sp <- ggscatter(df, x = "nExp", y = y_axis, color = "TumorType", palette = cols,
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray", linetype = 2), # Customize reg. line
                  conf.int = TRUE # Add confidence interval
  )
  print(sp + stat_cor(method = "spearman") + ggtitle("Spearman correlation") + labs(caption = select))
  
  
  sp <- ggscatter(df, x = "nExp", y = y, color = "TumorType", palette = cols,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )
  print(sp + stat_cor(aes(color = TumorType), method = "spearman") + ggtitle("Spearman correlation") + labs(caption = select))
dev.off()
