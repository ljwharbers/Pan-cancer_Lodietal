### Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

###### Accessing high-quality single-cell data
# To generate a Seurat object containing all shared pan-cancer high-quality cells, users can access the read count data for each individual cancer type on our lab’s website (https://lambrechtslab.sites.vib.be/en/dataaccess).
# After performing annotation analysis (refer to "Individual Cancer Types Analysis" and the cell subtype annotations in the "Cell Subtype Annotation" folder), users can explore the metadata of our pan-cancer atlas, available as "Lodietall_metadata.csv" within the "Master_files" folder. 
# This metadata file includes both general sample information and detailed cell subtype annotations for further reference.
# Once this object has been generated, proceed with the following analysis. 

object <- readRDS("/path/object.RDS")

df <- object@meta.data

### Create cell count cells
cell_count <- df %>%
  group_by(TumorType) %>%
  summarise(n=n())
ordered_tumor <- cell_count$TumorType[order(cell_count$n)]
df$TumorType <- factor(df$TumorType, levels = ordered_tumor)

### Get the cell counts
cell_count <- df %>%
  group_by(TumorType) %>%
  summarise(n=n()) %>%
  #  mutate(n_lab = case_when(n < 100000 ~ n + 25000,
  #                           n >= 100000 ~ n - 28000)) %>% 
  #  mutate(n_loc = case_when(n < 30000 ~ n + 15000,
  #                           n >= 100000 ~ n - 50000,
  #                           n >= 30000 ~ n - 18000)) %>%  
  mutate(n_loc = n + 18000) %>%
  mutate(n_lab = paste0(round(n/1000), "k")) %>% 
  ggplot(aes(x=n, y=TumorType)) + 
  geom_bar(stat="identity", color="black", fill="orange") + 
  geom_text(aes(x=n_loc , y=TumorType, label=n_lab)) + 
  labs(x = "Cell count", y=NULL) + 
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE), n.breaks = 5) + 
  theme(axis.text.y = element_text(size=14))

sample_count <- df %>%
  dplyr::select(TumorType, BiopsySite, SampleID) %>%
  unique() %>%
  group_by(TumorType, BiopsySite) %>%
  summarise(n=n()) %>%
  mutate(BiopsySite = factor(BiopsySite, 
                             levels = manuscript_levels$BiopsySite$levels, 
                             labels = manuscript_levels$BiopsySite$labels)) %>% 
  ggplot(aes(x=n, y=TumorType, fill=BiopsySite)) + 
  geom_bar(stat="identity", color="black") + 
  #geom_point(aes(x=n-3, y=TumorType), shape = 21, fill="white", size=9) + 
  geom_text(aes(y=TumorType, label=n), position=position_stack(vjust = 0.5)) + 
  labs(x = "Sample count", y = NULL) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="gray", size=0.1),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.4),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=14)) + NoLegend()

### Create cell fraction plot
cell_fraction <- df %>% 
  group_by(TumorType, BiopsySite) %>% 
  summarise(n=n()) %>%
  mutate(BiopsySite = factor(BiopsySite, 
                             levels = manuscript_levels$BiopsySite$levels, 
                             labels = manuscript_levels$BiopsySite$labels)) %>%
  ggplot(aes(x=n, y=TumorType, fill=BiopsySite)) +
  geom_bar(stat = "identity", position = "fill", color="black") + 
  labs(x = "Cell fraction", y = NULL, fill = "Biopsy site") + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(labels = scales::percent) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.4),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=14))


### Create figure of cell fractions and total cell counts
data = merge(dplyr::select(cell_count$data, -n), cell_fraction$data, by="TumorType")
cell_count = ggplot(data, aes(x=n, y=TumorType, fill=BiopsySite)) + 
  geom_bar(stat="identity", color="black") + 
  geom_text(aes(x=n_loc, y=TumorType, label=n_lab)) + 
  labs(x = "Cell count", y=NULL) + 
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE), n.breaks = 5) + 
  theme(axis.text.y = element_text(size=14)) + NoLegend() + xlim(0,180000) + 
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="gray", size=0.1),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.4),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=14)) + NoLegend()

### Create cell fraction plot
cell_types <- df %>% 
  group_by(TumorType, Majorcelltype_annotation) %>% 
  summarise(n=n()) %>%
  ggplot(aes(x=n, y=TumorType, fill=Majorcelltype_annotation)) +
  geom_bar(stat = "identity", position = "fill", color="black") + 
  labs(x = "Celltype fraction", y = NULL) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

### Create patient plot
patient_count <- df %>%
  dplyr::select(TumorType, Patient) %>%
  unique() %>%
  group_by(TumorType) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=n, y=TumorType)) + 
  geom_bar(stat="identity", color="black", fill="orange") + 
  geom_text(aes(x=n-3, y=TumorType, label=n)) + 
  labs(x = "Patient count", y = NULL) + 
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="gray", size=0.1),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.4),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=14))

### Pool the plots
p = ((cell_count + 
        sample_count) /
       (patient_count +
          cell_fraction)) + 
  guides(fill = guide_legend(override.aes=list(shape = 20)))
plot_layout(guides="collect") & theme(axis.text = element_text(color="black"))

### View plot
figsize(9, 12)
p

### Save plot
pdf("/path/Ncells_Nsamples_NPatients_Ncellfraction_acrosscancertypes.pdf", width = 9, height = 12)
p
dev.off()

