# ###### AddModules function to be added in shinyapp#######
# # readRDS("shinyApp/data.avg.rds") -> data.avg
# # readRDS("shinyApp/assay_data.rds") -> assay.data
# seu = AddModuleScore(seu, features = list(c("ITGAE", "PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1", "TNFRSF9", "TNFRSF18")), name="Tcell_Reactivity")
# seu$Tcell_Reactivity = seu$Tcell_Reactivity1
# seu$Tcell_Reactivity1 = NULL
# 
# 
# white2red <- c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
#                "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")
#   
# FeaturePlot(seu,features = "Tcell_Reactivity") +  scale_color_gradientn("", colors = white2red)

# nbin = 24
# ctrl = 100
# features <- c("ITGAE", "PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1", "TNFRSF9", "TNFRSF18")

# # Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations. The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
# data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
#                                 n = nbin,
#                                 labels = FALSE,
#                                 right = FALSE)
# # Set the names of the cuts as the gene names
# names(x = data.cut) <- names(x = data.avg)
# # Create an empty list the same length as the number of input gene sets. This will contain the names of the control genes
# ctrl.use <- vector(mode = "list", length = 1)
# for (j in 1:length(x = features)) {
#   ctrl.use[[1]] <- c(
#       ctrl.use[[1]],
#       names(x = sample(
#         x = data.cut[which(x = data.cut == data.cut[features[j]])],
#         size = ctrl,
#         replace = FALSE
#     ))
#     )
# }

# # Remove any repeated gene names - even though we set replace=FALSE when we sampled genes from the same expression bin, there may be more than two genes in our input gene list that fall in the same expression bin, so we can end up sampling the same gene more than once.
# ctrl.use <- lapply(X = ctrl.use, FUN = unique)

# ## Get control gene scores
# # Create an empty matrix with dimensions;
# # number of rows equal to the number of gene sets (just one here)
# # number of columns equal to number of cells in input Seurat object
# ctrl.scores <- matrix(
#   data = numeric(length = 1L),
#   nrow = 1,
#   ncol = ncol(assay.data)
# )
# features.use <- ctrl.use[[1]]
# ctrl.scores[1, ] <- Matrix::colMeans(x = assay.data[features.use, ])

# ## Get scores for input gene sets
# # Similar to the above, create an empty matrix
# features.scores <- matrix(
#   data = numeric(length = 1L),
#   nrow = 1,
#   ncol = ncol(x = assay.data)
# )

# data.use <- assay.data[features, , drop = FALSE]
# features.scores[1, ] <- Matrix::colMeans(x = data.use)

# features.scores.use <- features.scores - ctrl.scores

# rownames(x = features.scores.use) <- "Tcell_Reactivity"
# features.scores.use <- as.data.frame(x = t(x = features.scores.use))
# rownames(x = features.scores.use) <- colnames(x = assay.data)
# features.scores.use$sampleID <- row.names(features.scores.use)

# sc1meta = readRDS("shinyApp/sc1meta.rds")

# inner_join(features.scores.use,sc1meta,by="sampleID") -> sig1
# colnames(sig1)[2] <- "barcodes"
# sig1 <- left_join(sig1,barcodes)


######################### Calculcations dependent of score ###############################
# figsize = function(width=8, height=8) {
#   options(repr.plot.width=width, repr.plot.height=height)
# }
# figsize(9, 8)

# df2_major <- readRDS(file = "shinyApp/df2_major.rds")
# df_tumor <- readRDS(file = "shinyApp/df_tumor.rds")
# all_types <- readRDS(file = "shinyApp/alltypes.rds")
# counts <- readRDS(file = "shinyApp/counts.rds")

# df_react = sig1 %>%
#   group_by(SampleID) %>%
#   summarise_at(.vars = vars(Tcell_Reactivity), .funs = mean)

# df_react <- df_react[df_react$SampleID %in% row.names(df2_major),]
# # Merge the dataframes
# df2 = merge(df2_major, df_react, by.x = 0, by.y = "SampleID") %>% 
#   column_to_rownames("Row.names")


# df = merge(df2, counts, by=0) %>% column_to_rownames("Row.names")
# df$SampleID = rownames(df)
# df_copy = df
# df_orig = df

# ## I need to have cell types stored to use it
# cor_tcells = df[, c(all_types$tcell_types, "Tcell_Reactivity")] %>% 
#   cor_test(vars = c(all_types$tcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
#   filter(var2 == "Tcell_Reactivity") %>%
#   filter(var1 != "Tcell_Reactivity")

# # Keep only the samples with 10 or more B-cells
# df = filter(df_orig, `B cell_count` >= 10) 

# cor_bcells = df[, c(all_types$bcell_types, "Tcell_Reactivity")] %>% 
#   cor_test(vars = c(all_types$bcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
#   filter(var2 == "Tcell_Reactivity") %>%
#   filter(var1 != "Tcell_Reactivity")

# df = filter(df_orig, `Endothelial_count` >= 10)

# cor_endo = df[, c(all_types$endo_types[all_types$endo_types != "Specific"], "Tcell_Reactivity")] %>% 
#   cor_test(vars = c(all_types$endo_types[all_types$endo_types != "Specific"], "Tcell_Reactivity"), method = "spearman") %>% 
#   filter(var2 == "Tcell_Reactivity") %>%
#   filter(var1 != "Tcell_Reactivity")

# df = filter(df_orig, `Macrophages-Monocytes_count` >= 20)

# cor_myeloid = df[, c(all_types$macro_types, mono_types, "Neutrophils", "Tcell_Reactivity")] %>% 
#   cor_test(vars = c(all_types$macro_types, mono_types, "Neutrophils", "Tcell_Reactivity"), method = "spearman") %>% 
#   filter(var2 == "Tcell_Reactivity") %>%
#   filter(var1 != "Tcell_Reactivity")

# df = filter(df_orig, `DC_count` >= 10)

# cor_dc = df[, c(all_types$dc_types, "Tcell_Reactivity")] %>% 
#   cor_test(vars = c(all_types$dc_types, "Tcell_Reactivity"), method = "spearman") %>% 
#   filter(var2 == "Tcell_Reactivity") %>%
#   filter(var1 != "Tcell_Reactivity")

# ct_names = c(all_types$tcell_types, all_types$bcell_types, all_types$endo_types[all_types$endo_types != "Specific"], 
#              c(all_types$macro_types, all_types$mono_types, "Neutrophils"), all_types$dc_types)
# major_types = c(
#   rep("NK/T-cell", length(all_types$tcell_types)),
#   rep("B-cell", length(all_types$bcell_types)),
#   rep("Endothelial", length(all_types$endo_types[all_types$endo_types != "Specific"])),
#   rep("Macrophages/Monocytes", length(c(all_types$macro_types, all_types$mono_types, "Neutrophils"))),
#   rep("DC", length(all_types$dc_types))
# )

# # Create combined dataframe of all results
# rvalues = c(cor_tcells$cor, cor_bcells$cor, cor_endo$cor, cor_myeloid$cor, cor_dc$cor)
# pvalues = c(cor_tcells$p, cor_bcells$p, cor_endo$p, cor_myeloid$p, cor_dc$p)
# df_2 = data.frame(celltype = ct_names, Score = rvalues, Pval = pvalues, name = "Spearman", CellType_lev2.5 = major_types)

# figsize(30, 4)
# df_2 = add_significance(df_2, p.col = "Pval", output.col = "sig")
# df_2[df_2$sig == "ns", "sig"] = ""
# df_2[df_2$sig == "****", "sig"] = "***"

# pc_order = df_2 %>%
#   arrange(desc(Score)) %>% 
#   .$celltype
# df_2$celltype = factor(df_2$celltype, levels=pc_order)
# df_2$CellType_lev2.5 = factor(df_2$CellType_lev2.5, levels = c("NK/T-cell", "B-cell", "Macrophages/Monocytes", "DC", "Endothelial"))
# df_2$name = "Pan-cancer"

# # Draw the pancancer plot
# p_pancancer = ggplot(df_2) + 
#   geom_tile(aes(x = celltype, y = name, fill=Score), width=0.98, height=0.98) + 
#   geom_text(aes(x = celltype, y = name, label = sig), size=5) + 
#   facet_grid(~CellType_lev2.5, space = "free", scales = "free") + 
#   scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-1, 1)) + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         #axis.text.y = element_blank(),
#         panel.background = element_rect(fill="white", color="black"),
#         strip.text.x = element_text(size=14),
#         strip.background = element_rect(fill="white", color="black"),
#         axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5)) + 
#   labs(x = NULL, y = NULL)

# p_pancancer + theme(axis.text.y = element_text(size=16))

# pdf(paste0(figure_folder, "/PC2_reactivity_new_seu.pdf"), width=30, height=4)
# p_pancancer + theme(axis.text.y = element_text(size=16))
# dev.off()


# res = NULL

# ct_names =  c(all_types$tcell_types, all_types$bcell_types, all_types$endo_types[all_types$endo_types != "Specific"], 
#               c(all_types$macro_types, all_types$mono_types, "Neutrophils"), all_types$dc_types)
# major_types = c(
#   rep("NK/T-cell", length(all_types$tcell_types)),
#   rep("B-cell", length(all_types$bcell_types)),
#   rep("Endothelial", length(all_types$endo_types[all_types$endo_types != "Specific"])),
#   rep("Macrophages/Monocytes", length(c(all_types$macro_types, all_types$mono_types, "Neutrophils"))),
#   rep("DC", length(all_types$dc_types))
# )



# filter_samples <- TRUE
# # Iterate over the cancer types to calculate the correlations. Need to check if there are 
# # any samples that pass the filtering and create empty dataframe if there's not.
# for (cancer_type in unique(df_tumor$TumorType3)) {

#   rvalues = c()
#   pvalues = c()
#   df_orig = df_copy

#   # Keep only the samples of the selected tumor type
#   samples = df_tumor %>% 
#     filter(TumorType3 == cancer_type) %>%
#     .$SampleID %>%
#     as.vector()
#   df_orig = df_orig[samples, ]

#   # T-cells
#   df = if (filter_samples) filter(df_orig, `T cell_count` >= 20) else df_orig
#   if (nrow(df)) {
#     cor_tcells = df[, c(all_types$tcell_types, "Tcell_Reactivity")] %>% 
#       cor_test(vars = c(all_types$tcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
#       filter(var2 == "Tcell_Reactivity") %>%
#       filter(var1 != "Tcell_Reactivity")
#   } else {
#     cor_tcells = data.frame(cor = rep(NA, length(all_types$tcell_types)), 
#                             p = rep(NA, length(all_types$tcell_types)))
#   }

#   # B-cells
#   df = if (filter_samples) filter(df_orig, `B cell_count` >= 10) else df_orig
#   if (nrow(df)) {
#     cor_bcells = df[, c(all_types$bcell_types, "Tcell_Reactivity")] %>% 
#       cor_test(vars = c(all_types$bcell_types, "Tcell_Reactivity"), method = "spearman") %>% 
#       filter(var2 == "Tcell_Reactivity") %>%
#       filter(var1 != "Tcell_Reactivity")
#   } else {
#     cor_bcells = data.frame(cor = rep(NA, length(all_types$bcell_types)), 
#                             p = rep(NA, length(all_types$bcell_types)))
#   }

#   # Endothelial
#   df = if (filter_samples) filter(df_orig, `Endothelial_count` >= 10) else df_orig
#   if (nrow(df)) {
#     cor_endo = df[, c(all_types$endo_types[all_types$endo_types != "Specific"], "Tcell_Reactivity")] %>% 
#       cor_test(vars = c(all_types$endo_types[all_types$endo_types != "Specific"], "Tcell_Reactivity"), method = "spearman") %>% 
#       filter(var2 == "Tcell_Reactivity") %>%
#       filter(var1 != "Tcell_Reactivity")
#   } else {
#     cor_endo = data.frame(cor = rep(NA, length(all_types$endo_types[all_types$endo_types != "Specific"])), 
#                           p = rep(NA, length(all_types$endo_types[all_types$endo_types != "Specific"])))
#   }

#   # Myeloid
#   df = if (filter_samples) filter(df_orig, `Macrophages-Monocytes_count` >= 20) else df_orig
#   if (nrow(df)) {
#     cor_myeloid = df[, c(all_types$macro_types, all_types$mono_types, "Neutrophils", "Tcell_Reactivity")] %>% 
#       cor_test(vars = c(all_types$macro_types, all_types$mono_types, "Neutrophils", "Tcell_Reactivity"), method = "spearman") %>% 
#       filter(var2 == "Tcell_Reactivity") %>%
#       filter(var1 != "Tcell_Reactivity")
#   } else {
#     cor_myeloid = data.frame(cor = rep(NA, length(all_types$myeloid_types)), 
#                              p = rep(NA, length(all_types$myeloid_types)))
#   }

#   # DC
#   df = if (filter_samples) filter(df_orig, `DC_count` >= 10) else df_orig
#   if (nrow(df)) {
#     cor_dc = df[, c(all_types$dc_types, "Tcell_Reactivity")] %>% 
#       cor_test(vars = c(all_types$dc_types, "Tcell_Reactivity"), method = "spearman") %>% 
#       filter(var2 == "Tcell_Reactivity") %>%
#       filter(var1 != "Tcell_Reactivity")
#   } else {
#     cor_dc = data.frame(cor = rep(NA, length(all_types$dc_types)), 
#                         p = rep(NA, length(all_types$dc_types)))
#   }

#   # Create combined dataframe for the cancer type
#   rvalues = c(cor_tcells$cor, cor_bcells$cor, cor_endo$cor, cor_myeloid$cor, cor_dc$cor)
#   pvalues = c(cor_tcells$p, cor_bcells$p, cor_endo$p, cor_myeloid$p, cor_dc$p)
#   df_2 = data.frame(celltype = ct_names, Score = rvalues, Pval = pvalues, name = "Spearman", CellType_lev2.5 = major_types)
#   df_2 = add_significance(df_2, p.col = "Pval", output.col = "sig")
#   df_2[df_2$sig == "ns", "sig"] = ""
#   df_2[df_2$sig == "****", "sig"] = "***"
#   temp = df_2
#   temp$TumorType3 = cancer_type

#   # Merge with large dataframe
#   if(is.null(res))
#     res = temp
#   else 
#     res = rbind(res, temp)
# }

# figsize(30, 9)

# # Order by descending score of the pancancer data
# res$celltype = factor(res$celltype, levels=pc_order)
# res$CellType_lev2.5 = factor(res$CellType_lev2.5, levels = c("NK/T-cell", "B-cell", "Macrophages/Monocytes", "DC", "Endothelial"))
# res$TumorType3 = factor(res$TumorType3, levels = rev(c("BC_early", "BC_adv", "CC", "CRC", "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv")))

# # Draw the pancancer plot
# p = ggplot(res) + 
#   geom_tile(aes(x = celltype, y = TumorType3, fill = Score), width=0.98, height=0.98) + 
#   geom_text(aes(x = celltype, y = TumorType3, label = sig), size=5) + 
#   facet_grid(~CellType_lev2.5, space = "free", scales = "free") + 
#   scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-1, 1)) + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         #axis.text.y = element_blank(),
#         panel.background = element_rect(fill="white", color="black"),
#         strip.text.x = element_text(size=14),
#         strip.background = element_rect(fill="white", color="black"),
#         axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5)) + 
#   labs(x = NULL, y = NULL)

# p

# p_pancancer = p_pancancer + theme(axis.text.x = element_blank(),
#                                   axis.ticks.x = element_blank()) + NoLegend()

# figsize(30, 10)
# p_temp = p + theme(strip.background = element_blank(),
#                    strip.text.x = element_blank())
# p_temp = p_pancancer / p_temp + plot_layout(height=c(1, 7))
# p_temp = p_temp & theme(axis.text.y = element_text(size=16))

# figure_folder <- "figures_app/"
# pdf(paste0(figure_folder, "/PC2_reactivity_stratified_filtered_seurat.pdf"), width=30, height=10)
# p_temp
# dev.off()

# df = p$data %>% 
#   select(celltype, Score, TumorType3) %>%
#   pivot_wider(names_from = "celltype", values_from = "Score", values_fill = 0) %>%
#   column_to_rownames("TumorType3")

# pc_df = p_pancancer$data
# colnames(pc_df) = c("celltype", "Score", "Pval", "name", "CellType_lev2.5", "sig")
# pc_df$TumorType3 = "Pan-cancer"

# model <- hclust(dist(df))
# dhc <- as.dendrogram(model)
# # Rectangular lines

# ddata <- dendro_data(dhc, type = "rectangle")

# dp = ggplot(segment(ddata)) + 
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#   coord_flip() + 
#   scale_y_reverse(expand = c(0.2, 0)) + 
#   theme_void() +
#   theme(panel.background = element_blank(), 
#         plot.margin = margin(t = 20,  # Top margin
#                              r = 0,  # Right margin
#                              b = 40,  # Bottom margin
#                              l = 10))

# p$data$TumorType3 = factor(p$data$TumorType3, levels = model$labels[model$order])

# figsize(30, 10)
# p_pc = p_pancancer + 
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.x = element_blank(),
#         axis.text.y = element_text(size=16)) + 
#   NoLegend()

# p_pc$data$name = "Pan-cancer"
# p = p + theme(strip.text.x = element_blank(),
#               strip.background = element_blank(),
#               plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
#               axis.text.y = element_text(size=16))
# fp = plot_spacer() + p_pc + dp + p + plot_layout(widths = c(1, 10), heights = c(1, 10))
# fp = fp

# pdf(paste0(figure_folder, "/PC2_reactivity_stratified_dendro_filtered_seu.pdf"), width=30, height=10)
# fp
# dev.off()
