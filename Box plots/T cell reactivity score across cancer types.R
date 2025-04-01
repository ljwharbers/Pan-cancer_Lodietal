# Score T-cell Reactivity Signature Across Cancer Types
# Examine the expression of the T-cell reactivity signature at the sample level across different cancer types and rank the cancer types in ascending order based on signature values.

### Load libraries
library(Seurat)
library(dplyr)
library(ggpubr)

###### Accessing high-quality single-cell data
# To generate a Seurat object containing all shared pan-cancer high-quality cells, users can access the read count data for each individual cancer type on our lab’s website (https://lambrechtslab.sites.vib.be/en/dataaccess).
# After performing annotation analysis (refer to "Individual Cancer Types Analysis" and the cell subtype annotations in the "Cell Subtype Annotation" folder), users can explore the metadata of our pan-cancer atlas, available as "Lodietall_metadata.csv" within the "Master_files" folder. 
# This metadata file includes both general sample information and detailed cell subtype annotations for further reference.
# Once this object has been generated, proceed with the following analysis. 

object <- readRDS("/path/object.RDS")

### Define T-cell reactivity signature 
Tcell_Reactivity <- c("ITGAE", "PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1", "TNFRSF9", "TNFRSF18")
object <- AddModuleScore(object = object, features = list(Tcell_Reactivity),
                         ctrl = length(Tcell_Reactivity), name = 'Tcell_Reactivity')

### Get the metadata from the Seurat object
meta <- object@meta.data

### Keep only samples with more than 500 cells (see Methods)
meta <- meta %>% 
  group_by(SampleID) %>%
  mutate(n = n()) %>%
  filter(n >= 500)

### Further filter samples with >= 20 T-cells
samples <- meta %>%
  group_by(SampleID, Intermediatecelltype_annotation) %>%
  summarise(n = n(), .groups = 'drop') %>%
  filter(Intermediatecelltype_annotation == "T cell" & n >= 20) %>%
  pull(SampleID)

meta <- meta %>% filter(SampleID %in% samples)

### Define the score and features
score <- "Tcell_Reactivity1"
features <- c(score)
meta$Intermediatecelltype_annotation <- factor(meta$Intermediatecelltype_annotation)
meta$TumorType <- factor(meta$TumorType)

### Calculate patient mean score
res <- meta %>% 
  filter(Intermediatecelltype_annotation == "T cell") %>%
  group_by(SampleID) %>%
  summarise(across(all_of(features), list(mean = ~mean(., na.rm = TRUE)), .names = "mean_{col}"), .groups = 'drop')

### Add tumor types to the results
tt <- meta %>% select(SampleID, TumorType) %>% distinct()
res <- merge(tt, res, by = "SampleID")

### Calculate mean score per tumor type
mean_score <- res %>% 
  group_by(TumorType) %>%
  summarise(mean_score = mean(mean_Tcell_Reactivity1, na.rm = TRUE), .groups = 'drop')

### Set TumorType levels based on mean scores
res$TumorType <- factor(res$TumorType, levels = mean_score$TumorType[order(mean_score$mean_score)])

### Define color palette for tumor types
tumor_type_colors <- c(
  "BC_adv" = "#9966c7ff",
  "BC_early" = "#d8b7f9ff",
  "CC" = "#de5ea3ff",
  "CRC" = "#bd5b60ff",
  "GBM" = "#e49c6aff",
  "HCC" = "#8fc887ff",
  "HGSOC" = "#558cc6ff",
  "HNSCC" = "#f8c352ff",
  "MEL" = "#a0a1a2ff",
  "NSCLC_adv" = "#389b95ff",
  "NSCLC_early" = "#84dcd2ff"
)


### Plot with significance annotations
p2 <- res %>%
  ggboxplot(x = "TumorType", y = "mean_Tcell_Reactivity1", fill = "TumorType", outlier.shape = NA, add = "mean") +
  geom_jitter(size = 0.5, width = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black")) +
  labs(x = NULL,
       y = NULL,
       title = "Reactive T-cell Signature") +
  theme(plot.caption = element_text(hjust = 0.5, size = 12)) +
  scale_fill_manual(values = tumor_type_colors) +  # Add tumor type colors
  coord_flip()

### Perform the wilcox test
ttests <- compare_means(mean_Tcell_Reactivity1 ~ TumorType, data = res, method = "wilcox") %>%
  filter(p < 0.05)
#    filter(p.adj.signif != "ns")

extract_all_comparisons <- function(.tbl) {
  .tbl %>%
    dplyr::select(.data$group1, .data$group2) %>%
    purrr::transpose() %>%
    purrr::modify_depth(1, unlist)
}
p2 = p2 + stat_compare_means(method = "wilcox", 
                             step.increase = 0.02,
                             label = "p.signif", 
                             comparisons = extract_all_comparisons(ttests), 
                             tip.length = 0.005, 
                             vjust = 1.1)

### Save the plot
ggsave("/path/Boxplots_TcellReactivity_acrossCancertypes.pdf", plot = p2) 
