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

# For consistency, we created one file containing T-cell expansion data in early and HNSCC ("meta_object_expansion_only_T_M_noNA.csv") and another containing T-cell reactivity scores ("Tcell_Reactivity_SampleID.csv"), both calculated at the sample level.
# In addition, we also created a file with the metadata of our cells ("metaData_ncells_SampleID.csv), in particular containing the total number of cells per sample, as well as "ncells_T_NK_SampleID.csv", with the number of T- and NK-cells per sample.
# For more details, please refer to the files in the "Master_files" folder.

### Read data of the above-mentioned files 

# Read expansion data
meta_object_expansion <- read.csv(paste0(dir, "/Master_files/meta_object_expansion_only_T_M_noNA.csv")) 

# Read file with T-cell reactivity scores per sample
Tcell_reactivity <- read.csv(paste0(dir, "/Master_files/Tcell_Reactivity_SampleID.csv"))

# Read file with metadata 
df_ncells_meta <- read.csv(paste0(dir, "/Master_files/metaData_ncells_SampleID.csv"))
ncells_T_NK <- read.csv(paste0(dir, "/Master_files/ncells_T_NK_SampleID.csv"))
df_ncells_meta_T_NK <- dplyr::full_join(df_ncells_meta, ncells_T_NK)
df_ncells_meta_T_NK[is.na(df_ncells_meta_T_NK)] <- 0

### Combine data
reactivity_expansion <- dplyr::left_join(meta_object_expansion, Tcell_reactivity)
reactivity_expansion %>% dplyr::filter(is.na(avg_signature)) # 0

### Filter frequencies based on annotation and samples (samples with >= 500 samples and only PrimaryTumor&Metastasis)
sample_list <- readRDS(paste0(dir,"SVM_sample_lists_PC2_PC3.rds"))
samples <- sample_list[["PC2"]]["react_samples"]
reactivity_expansion %>% dplyr::filter(!SampleID %in% samples$react_samples) # 3 

samples <- df_ncells_meta %>% filter(ncells_noLQ_noD >= 500, ncells_noLQ_noD_obj_TC >= 20, BiopsySite %in% c("Tumor", "Metastasis")) %>% dplyr::pull(SampleID)
reactivity_expansion %>% dplyr::filter(!SampleID %in% samples) # 3

samples <- df_ncells_meta_T_NK %>% filter(ncells_noLQ_noD >= 500, ncells_noLQ_noD_obj_TC >= 20, BiopsySite %in% c("Tumor", "Metastasis")) %>% dplyr::pull(SampleID)
reactivity_expansion %>% dplyr::filter(!SampleID %in% samples) # 3

reactivity_expansion.filtered <- reactivity_expansion %>% dplyr::filter(SampleID %in% samples)

included <- "500_T_M_20TnoNK"


### Create and save the plot
df <- reactivity_expansion.filtered
df$expansion <- factor(df$expansion, levels = c("NE", "E"))

colors <- c("midnightblue", "red2")

pdf(paste0(out_dir, "/boxplot_EvsNE_reactivity_", included, ".pdf"), width = 3.5, height = 3.5)
plot <- ggplot(df, aes_string(x= "expansion", y= "avg_signature", color = "expansion")) + 
  scale_y_continuous(expand = expansion(mult = c(.01, .1))) +
  geom_boxplot(position=position_dodge(1), fill = NA) + 
  geom_jitter(shape=16,  position=position_dodge(1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  scale_color_manual(name = "expansion", values = colors) +
  theme(text=element_text(size=16)) + 
  stat_compare_means(label = "p.format", method = "wilcox.test", paired = F, show.legend = F, size = 4)
print(plot + labs(caption = included))
dev.off()

# Per cancer type
pdf(paste0(out_dir, "/boxplot_EvsNE_TumorType_reactivity_", included, ".pdf"), width = 5, height = 3.5)
plot <- ggplot(df, aes_string(x= "expansion", y= "avg_signature", color = "expansion")) + 
  scale_y_continuous(expand = expansion(mult = c(.01, .1))) +
  geom_boxplot(position=position_dodge(1), fill = NA) + 
  geom_jitter(shape=16,  position=position_dodge(1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  scale_color_manual(name = "expansion", values = colors) +
  theme(text=element_text(size=16)) + 
  facet_wrap(~TumorType) +
  stat_compare_means(label = "p.format", method = "wilcox.test", paired = F, show.legend = F, size = 4)
print(plot + labs(caption = included))
dev.off()


### Scatterplot 
cols <- c(
  HNSCC = "#F2B342",
  BC_early = "#E0CCEE"
)

pdf(paste0(out_dir, "correlation_nExp2_reactivity_", included, ".pdf"), width = 5, height = 5)

sp <- ggscatter(df, x = "nExp", y = "avg_signature", color = "TumorType", palette = cols,
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray", linetype = 2), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
print(sp + stat_cor(method = "spearman") + ggtitle("Spearman correlation") + labs(caption = included))
  
sp <- ggscatter(df, x = "nExp", y = "avg_signature", color = "TumorType", palette = cols,
                add = "reg.line",  # Add regressin line
                conf.int = TRUE # Add confidence interval
)
print(sp + stat_cor(aes(color = TumorType), method = "spearman") + ggtitle("Spearman correlation") + labs(caption = included))

dev.off()





