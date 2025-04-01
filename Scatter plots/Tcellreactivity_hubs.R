# Load libraries 
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Specify directories
dir <- "/path/"
out_dir <- "/path/"

# For consistency, we created one file containing cell proportions ("Freq_SampleID.csv") and another containing T-cell reactivity scores ("Tcell_Reactivity_SampleID.csv"), both calculated at the sample level.
# Specifically, in the first file, cell proportions are computed at the sample level:
# * For major cell types, proportions are calculated relative to the total number of cells.
# * At the cell subtype level, proportions are calculated relative to the sum of cells within the corresponding major cell type.
# In addition, we also have a file with the metadata of our cells ("metaData_ncells_SampleID.csv), in particular containing the total number of cells per sample. For more details, please refer to the files in the "Master_files" folder.

### Read data of the above-mentioned files 
# Read file with cell proportions
freq_all_long <- read.csv(paste0(dir, "/Master_files/freq_all_SampleID_long.csv"))

# Read file with metadata 
df_ncells_meta <- read.csv(paste0(dir, "Master_files/metaData_ncells_SampleID.csv"))

# Read file with T-cell reactivity scores per sample
Tcell_reactivity <- read.csv(paste0(dir, "/Master_files/Tcell_Reactivity_SampleID.csv"))

### Filter frequencies based on annotation and samples (samples with >= 500 samples and only PrimaryTumor&Metastasis)
level = "Minorcelltype_annotation"
freq.long <- freq_all_long %>% filter(lev == level)

samples <- df_ncells_meta %>% filter(ncells_noLQ_noD >= 500, BiopsySite %in% c("Tumor", "Metastasis")) %>% dplyr::pull(SampleID) %>% unique()
included <- "500_T_M"

freq.long <- freq.long %>% dplyr::filter(SampleID %in% samples)

#### Select cell subtypes belonging to Type-1 immunity hub
freq.long <- freq.long %>% dplyr::filter(CellType %in%
                                           c("CD4+ TREG", "mRegDC", "Mono-like Mac", "LAM2", "Neutrophils", "AXL_DC", "Inflam Mac", "Lymphatic",
                                             "CD4+ TH1", "CD8+ TEX", "Prolif T"))
select <- "Type-1 immunity hub"

### Select cell subtypes belonging to TLS-like hub
# freq.long <- freq.long %>% dplyr::filter(CellType %in%
#                                            c("CD4+ TFH", "GC B","Breg", "Plasmablast", "IFN Mac",
#                                              "IgG mature", "IgG immature", "IgA mature", "IgA immature", "mQuiescDC"))
# 
# 
# select <- "TLS-like hub"

### Add missing values
freq.long0 <- freq.long %>% tidyr::complete(SampleID, CellType, lev)
freq.long0[is.na(freq.long0)] <- 0

### Combine cell subtypes of the corresponding hub
freq_sum <- freq.long0 %>% dplyr::group_by(SampleID) %>% dplyr::summarise(n_sum = sum(n), freq_sum = sum(freq))

ncells <- df_ncells_meta %>% dplyr::select(SampleID, TumorType, ncells, ncells_noLQ_noD)

freq_sum <- dplyr::left_join(freq_sum, ncells)
freq_sum$combined_prop.clean <- freq_sum$n_sum / freq_sum$ncells_noLQ_noD

freq_sum <- dplyr::left_join(freq_sum, Tcell_reactivity)

### Select only samples with >= 500 samples and only PrimaryTumor&Metastasis. 
# For consistency, we created one file containing the number of T-cells and NK-cells ("ncells_T_NK_SampleID.csv"). For more details, please refer to the files in the "Master_files" folder.
ncells_T_NK <- read.csv(paste0(dir, "/path/ncells_T_NK_SampleID.csv"))
df_ncells_meta_T_NK <- dplyr::full_join(df_ncells_meta, ncells_T_NK)
df_ncells_meta_T_NK[is.na(df_ncells_meta_T_NK)] <- 0

samples <- df_ncells_meta_T_NK %>% filter(ncells_noLQ_noD >= 500, ncells_TC_noNK >= 20, BiopsySite %in% c("Tumor", "Metastasis"))%>% dplyr::pull(SampleID)

freq_sum_filtered <- freq_sum %>% dplyr::filter(SampleID %in% samples)
included <- "500_T_M_20TnoNK"

### Create and save the plot
df <- freq_sum_filtered

cols <- c(
  HGSOC = "#317EC2",
  CC = "#C43E96",
  CRC = "#C03830",
  GBM = "#E7872B",
  HNSCC = "#F2B342",
  HCC = "#94C47D",
  MEL = "#BBBBBC",
  BC_early = "#E0CCEE",
  BC_adv = "#9E66CD",
  NSCLC_early = "#BFE2E1",
  NSCLC_adv = "#46A19C")

#### Across cancer type
y <- c("combined_prop.clean")
pdf(paste0(out_dir, "Scatterplot_correlation_TcellReactivity_vs_", select, "_", included, ".pdf"), width = 6, height = 5)

sp <- ggscatter(df, x = "avg_signature", y = y, color = "TumorType3", palette = cols,
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray", linetype = 2), # Customize reg. line
                  conf.int = TRUE # Add confidence interval
  ) 
print(sp + stat_cor(method = "spearman") + ggtitle("Spearman correlation") + labs(caption = included) + theme(legend.position = "right"))
dev.off()

