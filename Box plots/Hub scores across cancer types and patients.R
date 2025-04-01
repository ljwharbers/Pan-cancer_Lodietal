### Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(forcats)
library(tibble)

# For consistency, we created one file containing cell proportions ("freq_all_SampleID_long.csv"). Specifically, in the this file, cell proportions are computed at the sample level:
# * For major cell types, proportions are calculated relative to the total number of cells.
# * At the cell subtype level, proportions are calculated relative to the sum of cells within the corresponding major cell type.
# In addition, we created a file with the metadata of our cells ("metaData_ncells_SampleID.csv), in particular containing the total number of cells per sample. For more details, please refer to the files in the "Master files" folder.


### Specify directories
dir <- "/path/"
out_dir <- "/path/"

### Read data of the above-mentioned files with cell proportion
# read file with cell proportions
freq_all_long <- read.csv(paste0(dir, "/freq_all_SampleID_long.csv"))

### Read file with metadata 
df_ncells_meta <- read.csv(paste0(dir, "/metaData_ncells_SampleID.csv"))

level = "Minorcelltype_annotation"
freq.long <- freq_all_long %>% filter(lev == level)

#### Select cell sucblusters belonging to Type-1 immunity hub
freq.long <- freq.long %>% dplyr::filter(CellType %in%
                                    c("CD4+ TREG", "mRegDC", "Mono-like Mac", "LAM2", "Neutrophils", "AXL_DC", "Inflam Mac", "Lymphatic",
                                      "CD4+ TH1", "CD8+ TEX", "Prolif T"))
hub <- "Type-1 immunity hub"

# #### Select cell sucblusters belonging to TLS-like hub
# freq.long <- freq.long %>% dplyr::filter(CellType %in%
#                                            c("CD4+ TFH", "GC B","Breg", "Plasmablast", "IFN Mac",
#                                              "IgG mature", "IgG immature", "IgA mature", "IgA immature", "mQuiescDC"))
# hub <- "TLS-like hub"


### Add missing values
freq.long0 <- freq.long %>% tidyr::complete(SampleID, CellType, lev)
freq.long0[is.na(freq.long0)] <- 0

### Combine cell subtypes of the corresponding hub
sample_freq_sum <- freq.long0 %>% dplyr::group_by(SampleID) %>% dplyr::summarise(n_sum = sum(n), freq_sum = sum(freq))

ncells <- df_ncells_meta %>% dplyr::select(SampleID, ncells, ncells_noLQ_noD, TumorType, BiopsySite)
sample_freq_sum <- dplyr::left_join(sample_freq_sum, ncells)
sample_freq_sum$combined_prop.clean <- sample_freq_sum$n_sum / sample_freq_sum$ncells_noLQ_noD

### Select only samples with >= 500 samples and only PrimaryTumor&Metastasis 
included <- "_500_PrimaryTumor&Metastasis"
select <- paste0(hub, "_500_PrimaryTumor&Metastasis")

df <- sample_freq_sum %>% filter(ncells_noLQ_noD >= 500, BiopsySite %in% c("Tumor", "Metastasis"))

### Create and save the plot 
colors <- c(
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
pdf(paste0(out_dir, "Boxplot_", select, "_TumorType_ranked.pdf"), width = 6, height = 5)

### Order by increasing mean per tumor type
mean_score <- df %>% dplyr::group_by(TumorType) %>% dplyr::summarise_at(vars(y), funs(mean(., na.rm=T)))
df$TumorType <- factor(df$TumorType, levels=mean_score$TumorType[order(mean_score[[y]])])

p2 <- df %>%
  dplyr::group_by(TumorType) %>%
  ggboxplot(x = "TumorType", y = y, fill="TumorType", outlier.shape = NA, add="mean")

p2 = p2 +
  geom_jitter(size=.5, width=0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=12, color = "black")) +
  theme(legend.position = "none",
        plot.title = element_text(size=20, face="bold", hjust=0.5),
        axis.text.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(angle=0, hjust=0.5, color = "black")) +
  labs(x=NULL,
       #y=NULL,
       title = hub) +
  theme(plot.caption = element_text(hjust = 0.5, size=12)) +
  scale_fill_manual(values = colors)  # Add tumor type colors

p2 = p2 + coord_flip()

### Perform the wilcox test
# ttests <- compare_means(combined_prop.clean ~ TumorType, data = df, method = "wilcox") %>%
#   filter(p < 0.05)
# 
# extract_all_comparisons <- function(.tbl) {
#   .tbl %>%
#     dplyr::select(.data$group1, .data$group2) %>%
#     purrr::transpose() %>%
#     purrr::modify_depth(1, unlist)
# }
# 
# p2 = p2 + stat_compare_means(method = "wilcox",
#                              step.increase = 0.02,
#                              label = "p.signif",
#                              comparisons = extract_all_comparisons(ttests),
#                              tip.length = 0.005,
#                              vjust = 1.1)

print(p2 + labs(caption = select))

dev.off()


### Across patients 
pdf(paste0(out_dir, "Boxplot", select, "_Patient_ranked.pdf"), width = 20, height = 5)

### Define the score to plot
score <- "combined_prop.clean"

### Order the data based on the selected score
df$SampleID <- factor(df$SampleID, levels = df[order(df[[score]]), ][["SampleID"]])

### Create the bar plot
p <- ggplot(data = df, aes(x = SampleID, y = !!sym(score), fill = TumorType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +  # Custom colors
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.background = element_blank()) +
  labs(x = NULL, title = select)

### Show and save the plot
print(p)
dev.off()
