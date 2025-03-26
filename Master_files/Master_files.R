# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read obj with all shared cells across cancer types
object <- readRDS("/path/object_allcells.rds") 

# Define directory
out_dir <- "/path/"

# Remove low quality (LD), Doublets, tissue-specific major cell types and subtypes
Idents(object) <- "Intermediatecelltype_annotation"
object_clean_cellsubtypes <- subset(object, idents = c("Doublets", "Low quality","Enteric glia", "Erythroblast", "Erythrocyte", "Muscle cell","Oligodendrocytes", "Sinusoidal"), invert = T)

object_clean_cellsubtypes$Intermediatecelltype_annotation_obj <- as.character(object_clean_cellsubtypes$Intermediatecelltype_annotation)
object_clean_cellsubtypes$Intermediatecelltype_annotation_obj <- ifelse(object_clean_cellsubtypes$Intermediatecelltype_annotation_obj %in% c("Tcell", "NK"), 
                                                               "NK/T cell", object_clean_cellsubtypes$Intermediatecelltype_annotation_obj)

obj_minor_annotation <- object_clean_cellsubtypes@meta.data %>% dplyr::select(Intermediatecelltype_annotation_obj, Minorcelltype_annotation) %>% unique()
write.csv(obj_minor_annotation, paste0(out_dir, "Minor_Intermediatecelltype_annotation.csv"), row.names = F)

obj_minor_annotation_2.0 <- object_clean_cellsubtypes@meta.data %>% dplyr::select(Intermediatecelltype_annotation_obj, Minorcelltype_annotation, Minorcelltype_annotation_2.0) %>% unique()
write.csv(obj_minor_annotation_2.0, paste0(out_dir, "Minor_2.0_Intermediatecelltype_annotation.csv"), row.names = F)

# # Calculate number of cells per SampleID
# ncells <- object@meta.data %>% dplyr::group_by(SampleID) %>% dplyr::summarise(ncells = n())
# ncells_noLQnoD <- object_clean_celltypes@meta.data %>% dplyr::group_by(SampleID) %>% dplyr::summarise(ncells_noLQnoD = n())
# 
# ncells_no_specific <- object_clean_cellsubtypes@meta.data %>% dplyr::filter(Intermediatecelltype_annotation_obj %in% c("TC", "BC", "MM", "DC3", "EC")) %>%
#   dplyr::group_by(SampleID, Intermediatecelltype_annotation_obj) %>% dplyr::summarise(ncells_noLQnoD = n())
# 
# object_clean_cellsubtypes@meta.data %>% dplyr::filter(Intermediatecelltype_annotation_obj %in% c("TC", "BC", "MM", "DC3", "EC")) %>% dim()
# 
# # ncells_no_specific %>% tidyr::pivot_wider(names_from = "Intermediatecelltype_annotation_obj", values_from = "ncells_noLQnoD") %>% head()
# 
# ncells_no_specific.wide <- ncells_no_specific %>% tidyr::pivot_wider(names_from = "Intermediatecelltype_annotation_obj", values_from = "ncells_noLQnoD") 
# new_colnames <- colnames(ncells_no_specific.wide)
# new_colnames[-1] <- paste0("ncells_noLQ_noD_obj_", new_colnames[-1])
# colnames(ncells_no_specific.wide) <- new_colnames
# 
# ncells_no_specific.wide[is.na(ncells_no_specific.wide)] <- 0
# 
# df_ncells <- dplyr::full_join(ncells, ncells_noLQnoD)
# df_ncells <- dplyr::full_join(df_ncells, ncells_no_specific.wide)
# df_ncells[is.na(df_ncells)] <- 0
# 
# write.csv(df_ncells, paste0(out_dir, "ncells_summary_SampleID.csv"), row.names = F)
# 
# df_ncells1 <- read.csv(paste0(out_dir, "/master_files_01022024/ncells_summary_SampleID.csv"))
# df_ncells <- read.csv(paste0(out_dir, "ncells_summary_SampleID.csv"))
# 
# stopifnot(all.equal(df_ncells, df_ncells1))
# df_ncells == df_ncells1
# 
# table(df_ncells == df_ncells1, useNA = 'ifany')
# identical(df_ncells,df_ncells1) # TRUE

# Create meta file
meta_samples <- object@meta.data %>% dplyr::select(SampleID, Patient, Technology, TumorType, BiopsySite) %>% unique() 

df_ncells <- df_ncells %>% dplyr::rename(ncells_noLQ_noD = ncells_noLQnoD)
df_ncells_meta <- dplyr::full_join(meta_samples, df_ncells)

write.csv(df_ncells_meta, paste0(out_dir, "metaData_ncells_SampleID.csv"), row.names = F)


##### Calculate cell (sub)type proportion ----
df <- object_clean_celltypes@meta.data

# Major cell type annotation
freq_major_long <- df %>%
  droplevels() %>%
  group_by(SampleID, Majorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n))  # 1068

freq_major <- freq_major_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Majorcelltype_annotation", values_from = "freq", values_fill = 0)

# Intermediate cell type annotation
freq_intermediate_long <- df %>%
  droplevels() %>%
  group_by(SampleID, Intermediatecelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n))  # 1602

freq_intermediate <- freq_intermediate_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Intermediatecelltype_annotation", values_from = "freq", values_fill = 0) # 178

# merge Majorcelltype_annotation and Intermediatecelltype_annotation
freq_major_long$lev <- "Majorcelltype_annotation"
freq_intermediate_long$lev <- "Intermediatecelltype_annotation"

freq_major_long <- freq_major_long %>% dplyr::rename(CellType = Majorcelltype_annotation)
freq_intermediate_long <- freq_intermediate_long %>% dplyr::rename(CellType = Intermediatecelltype_annotation)

freq_celltypes_long <- rbind(freq_major_long, freq_intermediate_long)
write.csv(freq_celltypes_long, paste0(out_dir, "Freq_SampleID.csv"), row.names = F) 



##### Calculate proportions cell subtypes ----
df <- object_clean_cellsubtypes@meta.data

# TCell - Minorcelltype_annotation
TC_freq_long <- df %>%
  filter(Intermediatecelltype_annotation_obj == "TCell") %>%
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) # 3150

TC_freq <- TC_freq_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0) # 175

# Bell - Minorcelltype_annotation
BC_freq_long <- df %>%
  filter(Intermediatecelltype_annotation_obj == "BCell") %>%
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) # 1620

BC_freq <- BC_freq_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0) # 162

# DC - Minorcelltype_annotation
DC_freq_long <- df %>%
  filter(Intermediatecelltype_annotation_obj == "DC") %>%
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) # 1384

DC_freq <- DC_freq_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0) # 173

# Macrophage/Monocyte - Minorcelltype_annotation
MM_freq_long <- df %>%
  filter(Intermediatecelltype_annotation_obj == "Macrophage/Monocyte") %>%
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) # 2655

MM_freq <- MM_freq_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0) # 177

# EC - Minorcelltype_annotation
EC_freq_long <- df %>%
  filter(Intermediatecelltype_annotation_obj == "EC") %>%
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) 

EC_freq <- EC_freq_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation", values_from = "freq", values_fill = 0)

# # TCell: Tex and Treg subtypes - Minorcelltype_annotation_2.0
# TC_ex_reg_freq_long <- df %>%
#   filter(Intermediatecelltype_annotation_obj == "TCell") %>%
#   droplevels() %>%
#   group_by(SampleID, Minorcelltype_annotation_2.0, .drop=F) %>%
#   summarise(n=n()) %>%
#   mutate(freq = n/sum(n)) # 4900
# 
# TC_ex_reg_freq <- TC_ex_reg_freq_long %>%
#   dplyr::select(-n) %>%
#   pivot_wider(names_from = "Minorcelltype_annotation_2.0", values_from = "freq", values_fill = 0) # 175
# 
# write.csv(TC_ex_reg_freq_long, paste0(out_dir, "freq_TC_ex_reg_SampleID_long.csv"), row.names = F)
# write.csv(TC_ex_reg_freq, paste0(out_dir, "freq_TC_ex_reg_SampleID.csv"), row.names = F)


### TEX-cell and TREG-cell subtypes - Minorcelltype_annotation_2.0
EX <- df %>% dplyr::filter(Minorcelltype_annotation == "CD8+ TEX") 
TREG <- df %>% dplyr::filter(Minorcelltype_annotation == "CD4+ TREG")

ncells_TEX <- object_clean_cellsubtypes@meta.data %>% dplyr::filter(Minorcelltype_annotation == "CD8+ TEX") %>% dplyr::group_by(SampleID) %>% dplyr::summarise(ncells_TEX = n()) 
ncells_TREG <- object_clean_cellsubtypes@meta.data %>% dplyr::filter(Minorcelltype_annotation == "CD4+ TREG") %>% dplyr::group_by(SampleID) %>% dplyr::summarise(ncells_TREG = n()) 

ncells_TexTreg <- dplyr::full_join(ncells_TEX, ncells_TREG)

ncells_TexTreg %>% filter(is.na(ncells_TREG) | is.na(ncells_TEX))
ncells_TexTreg[is.na(ncells_TexTreg)] <- 0

# write.csv(ncells_TexTreg, paste0(out_dir, "ncells_TexTreg_SampleID.csv"), row.names = F)

Tex_freq_long <- df %>%
  filter(Minorcelltype_annotation == "CD8+ TEX") %>%
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation_2.0, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) 

Tex_freq_long$Minorcelltype_annotation_2.0 <- paste0(Tex_freq_long$Minorcelltype_annotation_2.0, "_lev7_Tex") 

Tex_freq <- Tex_freq_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation_2.0", values_from = "freq", values_fill = 0) # 158


Treg_freq_long <- df %>%
  filter(Minorcelltype_annotation == "CD4+ TREG") %>%
  droplevels() %>%
  group_by(SampleID, Minorcelltype_annotation_2.0, .drop=F) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n)) 

Treg_freq_long$Minorcelltype_annotation_2.0 <- paste0(Treg_freq_long$Minorcelltype_annotation_2.0, "_lev7_Treg")

Treg_freq <- Treg_freq_long %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = "Minorcelltype_annotation_2.0", values_from = "freq", values_fill = 0) # 161

TexTreg_freq_long <- rbind(Tex_freq_long, Treg_freq_long) 
TexTreg_freq <- dplyr::full_join(Tex_freq, Treg_freq)

# write.csv(TexTreg_freq_long, paste0(out_dir, "freq_TexTreg_SampleID_long.csv"), row.names = F)
# write.csv(TexTreg_freq, paste0(out_dir, "freq_TexTreg_SampleID.csv"), row.names = F)

# merge cell subtypes
TC_freq_long$lev <- "Minorcelltype_annotation"
TC_freq_long <- TC_freq_long %>% dplyr::rename(CellType = Minorcelltype_annotation)

BC_freq_long$lev <- "Minorcelltype_annotation"
BC_freq_long <- BC_freq_long %>% dplyr::rename(CellType = Minorcelltype_annotation)

DC3_freq_long$lev <- "Minorcelltype_annotation"
DC3_freq_long <- DC3_freq_long %>% dplyr::rename(CellType = Minorcelltype_annotation)

MM_freq_long$lev <- "Minorcelltype_annotation"
MM_freq_long <- MM_freq_long %>% dplyr::rename(CellType = Minorcelltype_annotation)

EC_freq_long$lev <- "Minorcelltype_annotation"
EC_freq_long <- EC_freq_long %>% dplyr::rename(CellType = Minorcelltype_annotation)


TC_ex_reg_freq_long$lev <- "Minorcelltype_annotation_2.0_TC"
TC_ex_reg_freq_long <- TC_ex_reg_freq_long %>% dplyr::rename(CellType = Minorcelltype_annotation_2.0)

TexTreg_freq_long$lev <- "Minorcelltype_annotation_2.0_TexTreg"
TexTreg_freq_long <- TexTreg_freq_long %>% dplyr::rename(CellType = Minorcelltype_annotation_2.0)


freq_cellSUBtypes_long <- rbind(TC_freq_long, BC_freq_long) 
freq_cellSUBtypes_long <- rbind(freq_cellSUBtypes_long, DC3_freq_long) 
freq_cellSUBtypes_long <- rbind(freq_cellSUBtypes_long, MM_freq_long) 
freq_cellSUBtypes_long <- rbind(freq_cellSUBtypes_long, EC_freq_long) 

# write.csv(freq_cellSUBtypes_long, paste0(out_dir, "freq_cellSUBtypes_SampleID_long.csv"), row.names = F)

freq_cellSUBtypes_all_long <- rbind(freq_cellSUBtypes_long, TC_ex_reg_freq_long) 
freq_cellSUBtypes_all_long <- rbind(freq_cellSUBtypes_all_long, TexTreg_freq_long) 

# write.csv(freq_cellSUBtypes_all_long, paste0(out_dir, "freq_cellSUBtypes_all_SampleID_long.csv"), row.names = F) # latest 17/03/2024 (17108)


freq_cellSUBtypes <- dplyr::full_join(TC_freq, BC_freq) # 
freq_cellSUBtypes <- dplyr::full_join(freq_cellSUBtypes, DC3_freq) 
freq_cellSUBtypes <- dplyr::full_join(freq_cellSUBtypes, MM_freq) 
freq_cellSUBtypes <- dplyr::full_join(freq_cellSUBtypes, EC_freq) 

# write.csv(freq_cellSUBtypes, paste0(out_dir, "freq_cellSUBtypes_SampleID.csv"), row.names = F)

freq_cellSUBtypes_all <- dplyr::full_join(freq_cellSUBtypes, TC_ex_reg_freq) 
freq_cellSUBtypes_all <- dplyr::full_join(freq_cellSUBtypes_all, TexTreg_freq) 

# write.csv(freq_cellSUBtypes_all, paste0(out_dir, "freq_cellSUBtypes_all_SampleID.csv"), row.names = F)


# merge cell types and cell subtypes
freq_all_long <- rbind(freq_celltypes_long, freq_cellSUBtypes_all_long) 
# write.csv(freq_all_long, paste0(out_dir, "Freq_SampleID.csv"), row.names = F) # 19796 - 5  # latest 17/03/2024 (19778)

freq_all <- dplyr::full_join(freq_celltypes, freq_cellSUBtypes_all)
# write.csv(freq_all, paste0(out_dir, "Freq_SampleID_short.csv"), row.names = F) # 178 - 116



#### Calculate Tcell reactivity signature ----
# Read obj with all shared cells across cancer types
object <- readRDS("/path/object_allcells.rds") 

# Calculate Tcell reactivitysignature on all cells
signature <- list(c("ITGAE", "PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1", "TNFRSF9", "TNFRSF18"))
object <- AddModuleScore(object, features = signature, name = "reactive")

avg_sign <- object@meta.data %>% dplyr::filter(CellType_lev2.5 == "T cell") %>% 
  dplyr::select(SampleID, Tcell_Reactivity) %>% dplyr::group_by(SampleID) %>% dplyr::summarise(avg_signature = mean(Tcell_Reactivity)) # 175

write.csv(avg_sign, paste0(out_dir, "Tcell_Reactivity_SampleID.csv"), row.names = F)


### Calculate the number of T-cells, excluding NK-cells, per sample ----
# Read obj with all shared cells across cancer types
object <- readRDS("/path/object_allcells.rds") 

# Remove low quality (LD), Doublets, tissue-specific major cell types and subtypes
Idents(object) <- "Intermediatecelltype_annotation"
object_clean_cellsubtypes <- subset(object, idents = c("Doublets", "Low quality","Enteric glia", "Erythroblast", "Erythrocyte", "Muscle cell","Oligodendrocytes", "Sinusoidal"), invert = T)

ncells_TC <- object_clean_cellsubtypes@meta.data %>% dplyr::filter(Intermediatecelltype_annotation %in% c("T cell", "NK")) %>% dplyr::group_by(SampleID) %>% dplyr::summarise(ncells_TC = n())
ncells_TC_noNK <- object_clean_cellsubtypes@meta.data %>% dplyr::filter(Intermediatecelltype_annotation == "T cell") %>% dplyr::group_by(SampleID) %>% dplyr::summarise(ncells_TC_noNK = n())

ncells_T_NK <- dplyr::full_join(ncells_TC, ncells_TC_noNK)

# write.csv(ncells_T_NK, paste0(out_dir, "ncells_T_NK_SampleID.csv"), row.names = F)