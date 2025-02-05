library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(gridExtra)
library(dplyr)
library(stringr)
library(grid)
library(cowplot)
library(GSEABase)
library(shinycustomloader)
library(ggpubr)
library(ggsignif)
library(plotly)
library(ggforce)
library(zoo)
library(tidyverse)
library("Hmisc")
library(corrplot)
library(patchwork)
library(Seurat)
library(rstatix)
set.seed(1234)

sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
sc1gene = readRDS("sc1gene.rds")
sc1meta = readRDS("sc1meta.rds")

#Genes to plot gene signature
msigdbr_list <- readRDS("hallmarks.rds")

### Useful stuff 
# Colour palette 
cList = list(c("#0000FF", "#3838FF", "#7171FF", "#AAAAFF", "#E2E2FF", "#FFE2E2", "#FFAAAA", "#FF7171", "#FF3838", "#FF0000"),
              rev(c("#0000FF","#3F1FFC" ,"#4624FC", "#4C29FB" ,"#522DFA", "#5932FA", "#5F37F9",
                   "#653CF8", "#6B41F7", "#7147F6", "#774DF5" ,"#7E53F4" ,"#855AF2" ,"#8C62F0",
                   "#936BEE", "#9C76EC", "#A683E9", "#B193E5", "#BFABDF", "#D3D3D3")),
             c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154"),
             c("#000004FF", "#231151FF", "#5F187FFF", "#982D80FF", "#D3436EFF", "#F8765CFF", "#FEBA80FF", "#FCFDBFFF")) 
names(cList) = c("Blue-White-Red","Grey-Blue","White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple","Black-Violet-Yellow") 

# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "500px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(10,18,24) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 

# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  

# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 

figsize = function(width=8, height=8) {
  options(repr.plot.width=width, repr.plot.height=height)
}

# Correlation function
prepare_corrplot = function (mat, 
                             corr_metric = "spearman", 
                             plot_type = "full", 
                             diag = T, 
                             order = NULL, 
                             method = "color", 
                             addrect = NULL, 
                             insig = "blank", 
                             color = NULL, 
                             addgrid.col = NULL)
{
  res <- rcorr(as.matrix(mat), type = corr_metric)
  M <- res$r
  pmat <- res$P
  plot <- function() {
    corrplot(corr = M, p.mat = pmat, method = method, diag = diag, 
             type = plot_type, sig.level = c(0.001, 0.01, 0.05), 
             pch.cex = 0.9, insig = insig, pch.col = "black", 
             tl.col = "black", order = order, addrect = addrect, 
             col = color, addgrid.col = addgrid.col)
  }
  return(list(plot = plot, p = pmat, M = M))
}
figsize = function(width=8, height=8) {
  options(repr.plot.width=width, repr.plot.height=height)
}


### Common plotting functions 
## Function to plot gene signature

plot_gene_signature<- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                               inpsub1, inpsub2, inpH5, inpGene,
                               inpcols, inpfsz, inpcol,inpsiz,inplab, tab, save = FALSE){
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
  
  #Check if counts is already loaded. Save it glabally
  if(!exists("assay.data")){
    readRDS("data.avg.rds") ->> data.avg
    readRDS("assay_data.rds") ->> assay.data
    readRDS("df_prop.rds") ->> df_prop
    readRDS("subtypes.rds") ->> subtypes
    readRDS("major_types.rds") ->> major_types
    readRDS("SampleID_TumorType.rds") ->> SampleID_TumorType
    readRDS("barcode.rds") ->> barcodes
  }
  
  set.seed(1)
  features <- geneList$gene
  nbin = 24
  ctrl = 100
  # Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations. The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
  data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                  n = nbin,
                                  labels = FALSE,
                                  right = FALSE)
  # Set the names of the cuts as the gene names
  names(x = data.cut) <- names(x = data.avg)
  # Create an empty list the same length as the number of input gene sets. This will contain the names of the control genes
  ctrl.use <- vector(mode = "list", length = 1)
  for (j in 1:length(x = features)) {
    ctrl.use[[1]] <- c(
      ctrl.use[[1]],
      names(x = sample(
        x = data.cut[which(x = data.cut == data.cut[features[j]])],
        size = ctrl,
        replace = FALSE
      ))
    )
  }
  
  # Remove any repeated gene names - even though we set replace=FALSE when we sampled genes from the same expression bin, there may be more than two genes in our input gene list that fall in the same expression bin, so we can end up sampling the same gene more than once.
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  
  ## Get control gene scores
  # Create an empty matrix with dimensions;
  # number of rows equal to the number of gene sets (just one here)
  # number of columns equal to number of cells in input Seurat object
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = 1,
    ncol = ncol(assay.data)
  )
  features.use <- ctrl.use[[1]]
  ctrl.scores[1, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  
  ## Get scores for input gene sets
  # Similar to the above, create an empty matrix
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = 1,
    ncol = ncol(x = assay.data)
  )
  
  data.use <- assay.data[features, , drop = FALSE]
  features.scores[1, ] <- Matrix::colMeans(x = data.use)
  
  features.scores.use <- features.scores - ctrl.scores
  
  rownames(x = features.scores.use) <- "Gene_signature"
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = assay.data)
  features.scores.use$sampleID <- row.names(features.scores.use)
  
  inner_join(features.scores.use,sc1meta,by="sampleID") -> sig1
  colnames(sig1)[2] <- "barcodes"
  sig1 <- left_join(sig1,barcodes)
  
  
  #If subset button was used, filter rows that we are not interested, both for UMAP or further calculation
  bgCells = FALSE 
  sig1$plot <- TRUE
  if(!is.null(inpsub2)){
    bgCells = TRUE 
    sig1[,inpsub1] %in% inpsub2 -> sig1$plot
  }
  
  
  ## Use only the selected ones
  df_sign = sig1[sig1$plot == TRUE,] %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(Gene_signature = mean(Gene_signature), ncells = n())%>% dplyr::filter(ncells >= 20)
  
  # Just filter samples to make sure dataframes are the same size
  df_sign <- df_sign[df_sign$SampleID %in% df_prop$SampleID,]
  #Now generate final plot
  if(tab == "correlation"){
  df_orig <- dplyr::inner_join(df_sign, df_prop, by = "SampleID")
  df_copy = df_orig
  
  # CORRELATION FOR ALL M+T SAMPLES WITH MORE THAN 500 CELLS AND MORE THAN X (20 or 10) SIGNATURE CELLTYPE AND SUBTYPE SPECIFIC FILTERING 
  # Choose to filter samples
  filter_samples = TRUE
  
  # TC
  # Keep only the samples with 20 or more TCs  (object, so including NK)
  df = if (filter_samples) filter(df_orig, ncells_noLQ_noD_obj_TC >= 20) else df_orig
  
  cor_TC = df[, c(subtypes$TC_types, "Gene_signature")] %>% 
    cor_test(vars = c(subtypes$TC_types, "Gene_signature"), method = "spearman") %>% 
    filter(var2 == "Gene_signature") %>%
    filter(var1 != "Gene_signature")
  
  
  # BC
  # Keep only the samples with 10 or more BCs
  df = if (filter_samples) filter(df_orig, ncells_noLQ_noD_obj_BC >= 10) else df_orig
  
  cor_BC = df[, c(subtypes$BC_types, "Gene_signature")] %>% 
    cor_test(vars = c(subtypes$BC_types, "Gene_signature"), method = "spearman") %>% 
    filter(var2 == "Gene_signature") %>%
    filter(var1 != "Gene_signature")
  
  # MM
  # Keep only the samples with 20 or more MMs 
  df = if (filter_samples) filter(df_orig, ncells_noLQ_noD_obj_MM >= 20) else df_orig
  
  cor_MM = df[, c(subtypes$MM_types, "Gene_signature")] %>% 
    cor_test(vars = c(subtypes$MM_types, "Gene_signature"), method = "spearman") %>% 
    filter(var2 == "Gene_signature") %>%
    filter(var1 != "Gene_signature")
  
  
  # DC
  # Keep only the samples with 10 or more DCs
  df = if (filter_samples) filter(df_orig, ncells_noLQ_noD_obj_DC3 >= 10) else df_orig
  
  cor_DC3 = df[, c(subtypes$DC3_types, "Gene_signature")] %>% 
    cor_test(vars = c(subtypes$DC3_types, "Gene_signature"), method = "spearman") %>% 
    filter(var2 == "Gene_signature") %>%
    filter(var1 != "Gene_signature")
  
  
  # EC
  # Keep only the samples with 10 or more ECs
  df = if (filter_samples) filter(df_orig, ncells_noLQ_noD_obj_EC >= 10) else df_orig
  
  cor_EC = df[, c(subtypes$EC_types, "Gene_signature")] %>% 
    cor_test(vars = c(subtypes$EC_types, "Gene_signature"), method = "spearman") %>% 
    filter(var2 == "Gene_signature") %>%
    filter(var1 != "Gene_signature")
  
  
  # Create combined dataframe of all results
  cor_subtypes <- bind_rows(list(cor_TC, cor_BC, cor_MM, cor_DC3, cor_EC))
  cor_df = cor_subtypes %>% dplyr::rename(Score = cor, Pval = p, name = method, subtype = var1) %>% dplyr::left_join(major_types, by = "subtype") %>% dplyr::select(-var2)
  
  
  # PLOT PANCANCER
  cor_df = add_significance(cor_df, p.col = "Pval", output.col = "sig")
  cor_df[cor_df$sig == "ns", "sig"] = ""
  cor_df[cor_df$sig == "****", "sig"] = "***"
  
  pc_order = cor_df %>%
    arrange(desc(Score)) %>% 
    .$subtype
  cor_df$subtype = factor(cor_df$subtype, levels=pc_order)
  cor_df$major_type = factor(cor_df$major_type, levels = c("NK/T cells", "B cells",  "Monocytes and macrophages", "Dendritic cells", "Endothelial cells"))
  cor_df$name = "Pan-cancer"
  
  # Draw the pancancer plot
  figsize(30, 4)
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
  
  
  
  # STRATIFIED PER TUMORTYPE - CORRELATION FOR ALL M+T SAMPLES WITH MORE THAN 500 CELLS AND MORE THAN  X (20 or 10) SIGNATURE CELLTYPE AND SUBTYPE SPECIFIC FILTERING 
  res = NULL
  
  # Iterate over the cancer types to calculate the correlations. Need to check if there are 
  # any samples that pass the filtering and create empty dataframe if there's not.
  
  for(cancer_type in as.character(unique(sig1$Cancer_type_early_adv[sig1$plot == T]))) {
    
    # for (cancer_type in unique(SampleID_TumorType$TumorType3)) {
    
    
    # Keep only the samples of the selected tumor type
    samples = SampleID_TumorType %>% 
      dplyr::filter(TumorType3 == cancer_type) %>%
      .$SampleID %>%
      as.vector()
    
    df_orig_cancerType = df_orig %>% dplyr::filter(SampleID %in% samples)
    
    ###This line was added because if you have only one 1 sample the correlation gives an error
    if(nrow(df_orig_cancerType) > 2){
      # TC
      df = if (filter_samples) filter(df_orig_cancerType, ncells_noLQ_noD_obj_TC >= 20) else df_orig_cancerType
      if (nrow(df)) {
        cor_TC = df[, c(subtypes$TC_types, "Gene_signature")] %>% 
          cor_test(vars = c(subtypes$TC_types, "Gene_signature"), method = "spearman") %>% 
          filter(var2 == "Gene_signature") %>%
          filter(var1 != "Gene_signature")
      } else {
        cor_TC = data.frame(var1 = subtypes$TC_types,
                            var2 = "Gene_signature",
                            cor = rep(NA, length(subtypes$TC_types)),
                            statistic = rep(NA, length(subtypes$TC_types)),
                            p = rep(NA, length(subtypes$TC_types)),
                            method = "spearman")
      }
      
      
      
      # BC
      df = if (filter_samples) filter(df_orig_cancerType, ncells_noLQ_noD_obj_BC >= 10) else df_orig_cancerType
      if (nrow(df)) {
        cor_BC = df[, c(subtypes$BC_types, "Gene_signature")] %>% 
          cor_test(vars = c(subtypes$BC_types, "Gene_signature"), method = "spearman") %>% 
          filter(var2 == "Gene_signature") %>%
          filter(var1 != "Gene_signature")
      } else {
        cor_BC = data.frame(var1 = subtypes$BC_types,
                            var2 = "Gene_signature",
                            cor = rep(NA, length(subtypes$BC_types)),
                            statistic = rep(NA, length(subtypes$BC_types)),
                            p = rep(NA, length(subtypes$BC_types)),
                            method = "spearman")
      }
      
      
      # EC
      df = if (filter_samples) filter(df_orig_cancerType, ncells_noLQ_noD_obj_EC >= 10) else df_orig_cancerType
      if (nrow(df)) {
        cor_EC = df[, c(subtypes$EC_types, "Gene_signature")] %>% 
          cor_test(vars = c(subtypes$EC_types, "Gene_signature"), method = "spearman") %>% 
          filter(var2 == "Gene_signature") %>%
          filter(var1 != "Gene_signature")
      } else {
        cor_EC = data.frame(var1 = subtypes$EC_types,
                            var2 = "Gene_signature",
                            cor = rep(NA, length(subtypes$EC_types)),
                            statistic = rep(NA, length(subtypes$EC_types)),
                            p = rep(NA, length(subtypes$EC_types)),
                            method = "spearman")
      }
      
      
      
      # MM
      df = if (filter_samples) filter(df_orig_cancerType, ncells_noLQ_noD_obj_MM >= 20) else df_orig_cancerType
      if (nrow(df)) {
        cor_MM = df[, c(subtypes$MM_types, "Gene_signature")] %>% 
          cor_test(vars = c(subtypes$MM_types, "Gene_signature"), method = "spearman") %>% 
          filter(var2 == "Gene_signature") %>%
          filter(var1 != "Gene_signature")
      } else {
        cor_MM = data.frame(var1 = subtypes$MM_types,
                            var2 = "Gene_signature",
                            cor = rep(NA, length(subtypes$MM_types)),
                            statistic = rep(NA, length(subtypes$MM_types)),
                            p = rep(NA, length(subtypes$MM_types)),
                            method = "spearman")
      }
      
      
      
      
      # DC
      df = if (filter_samples) filter(df_orig_cancerType, ncells_noLQ_noD_obj_DC3 >= 10) else df_orig_cancerType
      if (nrow(df)) {
        cor_DC3 = df[, c(subtypes$DC3_types, "Gene_signature")] %>% 
          cor_test(vars = c(subtypes$DC3_types, "Gene_signature"), method = "spearman") %>% 
          filter(var2 == "Gene_signature") %>%
          filter(var1 != "Gene_signature")
      } else {
        cor_DC3 = data.frame(var1 = subtypes$DC3_types,
                             var2 = "Gene_signature",
                             cor = rep(NA, length(subtypes$DC3_types)),
                             statistic = rep(NA, length(subtypes$DC3_types)),
                             p = rep(NA, length(subtypes$DC3_types)),
                             method = "spearman")
      }
      
      
      
      # Create combined dataframe for the cancer typ
      cor_subtypes <- bind_rows(list(cor_TC, cor_BC, cor_MM, cor_DC3, cor_EC))
      cor_df = cor_subtypes %>% dplyr::rename(Score = cor, Pval = p, subtype = var1) %>% dplyr::left_join(major_types, by = "subtype") %>% dplyr::select(-var2)
      
      cor_df = add_significance(cor_df, p.col = "Pval", output.col = "sig")
      cor_df[cor_df$sig == "ns", "sig"] = ""
      cor_df[cor_df$sig == "****", "sig"] = "***"
      
      temp = cor_df
      temp$TumorType3 = cancer_type
      
      
      # Merge with large dataframe
      if(is.null(res))
        res = temp
      else 
        res = rbind(res, temp)
    }
    
  }
  
  # PLOT PANCANCER AND PER TUMORTYPE 
  # Order by descending score of the pancancer data
  res$subtype = factor(res$subtype, levels=pc_order)
  res$major_type = factor(res$major_type, levels = c("NK/T cells", "B cells",  "Monocytes and macrophages", "Dendritic cells", "Endothelial cells"))
  res$TumorType3 = factor(res$TumorType3, levels = rev(c("BC_early", "BC_adv", "CC", "CRC", "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv")))
  
  # Draw the pancancer plot
  figsize(30, 9)
  p = ggplot(res) + 
    geom_tile(aes(x = subtype, y = TumorType3, fill = Score), width=0.98, height=0.98) + 
    geom_text(aes(x = subtype, y = TumorType3, label = sig), size=5) + 
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
  
  
  # Add dendrogram
  df = p$data %>% 
    dplyr::select(subtype, Score, TumorType3) %>%
    pivot_wider(names_from = "subtype", values_from = "Score", values_fill = 0) %>%
    tibble::column_to_rownames("TumorType3")
  
  
  pc_df = p_pancancer$data
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
  
  # fp = plot_spacer() + p_pc + dp + p + plot_layout(widths = c(1, 10), heights = c(1, 10))
  
  #####    #Now join plots
  fp = wrap_plots(plot_spacer() + p_pc + dp + p + plot_layout(widths = c(1, 10), heights = c(1, 10)))
  last_generated <<- fp
  return(fp)
    
  } else if(tab =="expression"){
    ##################  Barplots ##################
    if(inpGrp == "Cancer_type_early_adv"){
      # Get colors from input
      ggCol = strsplit(inpConf[UI == inpGrp]$fCL, "\\|")[[1]]
      names(ggCol) = strsplit(inpConf[UI == inpGrp]$fID, "\\|")[[1]] 
      
      #Generate bar plot
      str_remove_all(df_sign$SampleID,"_[0-9]*$") %>% str_replace_all(.,"OV","HGSOC") -> df_sign$TumorType3
      
      temp <- df_sign %>%
        group_by(TumorType3) %>%
        summarise_at(.vars = vars(Gene_signature), .funs = mean) %>%
        arrange(Gene_signature) %>% .$TumorType3
      df_sign$TumorType3 <- factor(df_sign$TumorType3,levels=temp)
      
      barplot_p <- ggplot(df_sign,aes(x=Gene_signature,y=TumorType3,fill=TumorType3))+geom_boxplot()+geom_point()+
        scale_fill_manual(values = ggCol[temp]) +sctheme(base_size = 12, Xang = 45, XjusH = 1)+NoLegend()+ylab("Cancer type (early/adv)") + xlab("Gene signature (avg/sample)")
    }else{
      # Get colors from input
      ggCol = strsplit(inpConf[UI == inpGrp]$fCL, "\\|")[[1]]
      names(ggCol) = strsplit(inpConf[UI == inpGrp]$fID, "\\|")[[1]] 
      
      df_sign = sig1[sig1$plot == TRUE,] %>%
        group_by(SampleID,Cell_subtype) %>%
        summarise_at(.vars = vars(Gene_signature), .funs = mean)
      
      temp <- df_sign %>%
        group_by(Cell_subtype) %>%
        summarise_at(.vars = vars(Gene_signature), .funs = mean) %>%
        arrange(Gene_signature) %>% .$Cell_subtype
      df_sign$Cell_subtype <- factor(df_sign$Cell_subtype,levels=temp)
      
      barplot_p <- ggplot(df_sign,aes(x=Gene_signature,y=Cell_subtype,fill=Cell_subtype))+geom_boxplot()+geom_point()+
        scale_fill_manual(values = ggCol[temp]) +sctheme(base_size = 12, Xang = 45, XjusH = 1)+ NoLegend()+ylab("Cell subtype") + xlab("Gene signature (avg/sample/cell subtype)")
    }
    
    
    #####Generate UMAP
    ggOut = ggplot(sig1[sig1$plot == TRUE,], aes(UMAP_1, UMAP_2))
    if(bgCells){ 
      #Cells that were filtered out need to be shown in another color
      ggOut = ggOut + 
        geom_point(data = sig1[sig1$plot == FALSE,], color = "snow2", size = inpsiz,shape = 16)
    } 
    ggOut = ggOut + 
      geom_point(shape = 16, size = inpsiz, aes(color = Gene_signature)) + 
      xlab("UMAP1") + ylab("UMAP2") + 
      scale_color_gradientn("", colours = cList[[inpcol]])+
      guides(color = guide_colorbar(barwidth = 15))+ sctheme(base_size = sList[inpfsz])+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    
    if(inplab){ 
      ##Add group label # Change for every app
      label_df <- sig1 %>% group_by(Cell_subtype) %>% 
        summarise_at(.vars = vars(UMAP_1,UMAP_2),
                     .funs = c(mean="mean"))
      colnames(label_df) <- c("val","X","Y")
      lListX = min(nchar(paste0(label_df$val, collapse = "")), 200)
      lListX = lList - (0.25 * floor(lListX/50))
      ggOut = ggOut +
        geom_text_repel(data = label_df, aes(X, Y, label = val),
                        color = "grey10", bg.color = "grey95", bg.r = 0.15,
                        size = lListX[inpfsz], seed = 42)
    }
    
    
    #Now join plots
    blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
    geneset_plot <- grid.arrange(barplot_p,blank,ggOut,widths=c(0.55, 0.05, 0.4),nrow=1)
    last_generated_bar <<- geneset_plot
    return(geneset_plot)}
}




# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "val", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Do factoring if required 
  if(!is.na(inpConf[UI == inp1]$fCL)){ 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$val) 
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
    ggData$val = factor(ggData$val, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = sListX[inpfsz])) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
      lListX = lList - (0.25 * floor(lListX/50)) 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = lListX[inpfsz], seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 

scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
} 

# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  ### Modify to plot more than 1 per tab
  #Create a backup that will be used by every gene
  ggData -> ggData_backup
  shiny::validate(need(length(inp1) <= 8, "More than 8 genes to plot! Please reduce the gene list."))
  shiny::validate(need(length(inp1) > 0, "Please input at least 1 genes to plot"))
  #Generate 1 plot per gene
  for(i in 1:length(inp1)){
    #Load backup
    ggData <- ggData_backup
    #Gene of interest
    inp1[i] -> genes2plot
    #Load  H5
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    
    ggData$val = h5data$read(args = list(inpGene[genes2plot], quote(expr=))) 
    ggData[val < 0]$val = 0 
    h5file$close_all() 
    bgCells = FALSE 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
      bgCells = TRUE 
      ggData2 = ggData[!sub %in% inpsub2] 
      ggData = ggData[sub %in% inpsub2] 
    } 
    if(inpord == "Max-1st"){ 
      ggData = ggData[order(val)] 
    } else if(inpord == "Min-1st"){ 
      ggData = ggData[order(-val)] 
    } else if(inpord == "Random"){ 
      ggData = ggData[sample(nrow(ggData))] 
    } 
    
    # Actual ggplot 
    ggOut = ggplot(ggData, aes(X, Y, color = val)) 
    if(bgCells){ 
      ggOut = ggOut + 
        geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
    } 
    ggOut = ggOut + 
      geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
      sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
      scale_color_gradientn(genes2plot, colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
    if(inpasp == "Square") { 
      ggOut = ggOut + coord_fixed(ratio = rat) 
    } else if(inpasp == "Fixed") { 
      ggOut = ggOut + coord_fixed() 
    } 
    assign(paste0("p",i),ggOut)
  }
  
  if(i==2){
    ggOut = grid.arrange(p1,p2,ncol=2)
  }else if(i==3){
    ggOut = grid.arrange(p1,p2,p3,ncol=3)
  }else if(i==4){
    ggOut = grid.arrange(p1,p2,p3,p4,ncol=4)
  }
  else if(i==5){
    ggOut = grid.arrange(p1,p2,p3,p4,p5,ncol=4)
  }
  else if(i==6){
    ggOut = grid.arrange(p1,p2,p3,p4,p5,p6,ncol=4)
  }
  else if(i==7){
    ggOut = grid.arrange(p1,p2,p3,p4,p5,p6,p7,ncol=4)
  }
  else if(i==8){
    ggOut = grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=4)
  }
  
  return(ggOut) 
} 

# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Map colours 
  ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
  ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
  ggData$v0 = ggData$v1 + ggData$v2 
  ggData = gg[ggData, on = c("v1", "v2")] 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(v0)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-v0)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) + 
    scale_color_gradientn(inp1, colours = cList[[1]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 

scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = sList[inpfsz], XYval = TRUE) 
  return(ggOut) 
} 

scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpsub1, inpsub2, inpH5, inpGene){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 

# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inptyp, inppts, inpsiz, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  
  #### Make a loop to generate plots for more than one gene without having to modify the code too much
  shiny::validate(need(length(inp2) <= 6, "More than 6 genes to plot! Please reduce the gene list."))
  shiny::validate(need(length(inp2) > 0, "Please input at least 1 genes to plot"))
  
  for(i in 1:length(inp2)){
    inp2[i] -> genes2plot
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                     with = FALSE] 
    colnames(ggData) = c("X", "sub") 
    
    # Load in either cell meta or gene expr
    if(genes2plot %in% inpConf$UI){ 
      ggData$val = inpMeta[[inpConf[UI == genes2plot]$ID]] 
    } else { 
      h5file <- H5File$new(inpH5, mode = "r") 
      h5data <- h5file[["grp"]][["data"]] 
      ggData$val = h5data$read(args = list(inpGene[genes2plot], quote(expr=))) 
      ggData[val < 0]$val = 0 
      set.seed(42) 
      tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
      ggData$val = ggData$val + tmpNoise 
      h5file$close_all() 
    } 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
      ggData = ggData[sub %in% inpsub2] 
    } 
    
    # Do factoring 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$X) 
    ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
    ggData$X = factor(ggData$X, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
    
    # Actual ggplot 
    if(inptyp == "violin"){ 
      ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
    } else { 
      ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
    } 
    if(inppts){ 
      ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
    } 
    ggOut = ggOut + xlab(inp1) + ylab(genes2plot) + 
      sctheme(base_size = sList[inpfsz], Xang = 90, XjusH = 1) +  
      scale_fill_manual("", values = ggCol) +
      theme(legend.position = "none")
    assign(paste0("p",i),ggOut)
    
  }
  
  if(i==2){
    ggOut = grid.arrange(p1,p2,ncol=2)
  }else if(i==3){
    ggOut = grid.arrange(p1,p2,p3,ncol=2)
  }else if(i==4){
    ggOut = grid.arrange(p1,p2,p3,p4,ncol=2)
  }
  else if(i==5){
    ggOut = grid.arrange(p1,p2,p3,p4,p5,ncol=2)
  }
  else if(i==6){
    ggOut = grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
  }
  return(ggOut) 
} 

# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
  .SD[,.(pctCells = 100 * sum(nCells) / tot, 
         nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 

# Get gene list 
scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 

# Plot gene expression bubbleplot / heatmap 
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
  
  # Prepare ggData 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData = data.table() 
  for(iGene in geneList$gene){ 
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(tmp) = c("sampleID", "sub") 
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
    tmp$geneName = iGene 
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
    ggData = rbindlist(list(ggData, tmp)) 
  } 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
  
  # Aggregate 
  ggData$val = expm1(ggData$val) 
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)), 
                  by = c("geneName", "grpBy")] 
  ggData$val = log1p(ggData$val) 
  
  # Scale if required 
  colRange = range(ggData$val) 
  if(inpScl){ 
    ggData[, val:= scale(val), keyby = "geneName"] 
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  } 
  
  # hclust row/col if necessary 
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
  tmp = ggMat$geneName 
  ggMat = as.matrix(ggMat[, -1]) 
  rownames(ggMat) = tmp 
  if(inpRow){ 
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
    ggRow = ggplot() + coord_flip() + 
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                         labels = unique(ggData$grpBy), expand = c(0, 0)) + 
      scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                         labels = hcRow$labels$label, expand = c(0, 0.5)) + 
      sctheme(base_size = sList[inpfsz]) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
  } 
  if(inpCol){ 
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
    ggCol = ggplot() + 
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                         labels = hcCol$labels$label, expand = c(0.05, 0)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                         labels = unique(ggData$geneName), expand=c(0,0)) + 
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(), 
            axis.text.y = element_text(color = "white")) 
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
  } 
  
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){ 
    # Bubbleplot 
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
      geom_point() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_size_continuous("proportion", range = c(0, 8), 
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(color = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else { 
    # Heatmap 
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
      geom_tile() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(fill = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank()) 
  } 
  
  # Final tidy 
  ggLeg = g_legend(ggOut) 
  ggOut = ggOut + theme(legend.position = "none") 
  if(!save){ 
    if(inpRow & inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                   layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                   layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                   layout_matrix = rbind(c(1),c(2)))  
    }  
  } else { 
    if(inpRow & inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                  layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                  layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                  layout_matrix = rbind(c(1),c(2)))  
    }  
  } 
  return(ggOut) 
} 





### Start server code 
shinyServer(function(input, output, session) { 
  ### For all tags and Server-side selectize 
  observe_helpers() 
  optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE,
                       selected = "FSCN1", options = list(
                         maxOptions = 300,create = TRUE, persist = TRUE, render = I(optCrt)))
  updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE, 
                       selected = c("CD68","CLEC9A","CCR7","TCF4"), options = list( 
                         maxOptions = 300, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE, 
                       selected = "CLEC9A", options = list( 
                         maxOptions = 300, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE, 
                       selected = "CLEC9A", options = list( 
                         maxOptions = 300, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE, 
                       selected = "FSCN1", options = list( 
                         maxOptions = 300, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1c1inp2", server = TRUE, 
                       choices = names(sc1gene), 
                       selected = c("FSCN1","CLEC9A","CCR7","TCF4"), options = list( 
                         maxOptions = length(sc1conf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
  
  
  ### Plots for tab a1 
  output$sc1a1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a1oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,  
             input$sc1a1sub1, input$sc1a1sub2, 
             input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) 
  }) 
  output$sc1a1oup1.ui <- renderUI({ 
    withLoader(plotOutput("sc1a1oup1", height = pList[input$sc1a1psz]) ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1a1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
    }) 
  output$sc1a1oup1.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
    }) 
  output$sc1a1.dt <- renderDataTable({ 
    ggData = scDRnum(sc1conf, sc1meta, input$sc1a1inp1, input$sc1a1inp2, 
                     input$sc1a1sub1, input$sc1a1sub2, 
                     "sc1gexpr.h5", sc1gene, "Decile")  #input$sc1a1splt
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
  
  output$sc1a1oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
             input$sc1a1sub1, input$sc1a1sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) 
  }) 
  output$sc1a1oup2.ui <- renderUI({ 
    withLoader(plotOutput("sc1a1oup2", height = pList[input$sc1a1psz]),type = "html",loader = "dnaspin") 
  }) 
  output$sc1a1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
    }) 
  output$sc1a1oup2.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
    }) 
  
  
  ### Plots for tab a2 
  output$sc1a2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a2oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) 
  }) 
  output$sc1a2oup1.ui <- renderUI({ 
    withLoader(plotOutput("sc1a2oup1", height = pList[input$sc1a2psz]) ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1a2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
    }) 
  output$sc1a2oup1.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
    }) 
  
  output$sc1a2oup2 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) 
  }) 
  output$sc1a2oup2.ui <- renderUI({ 
    withLoader(plotOutput("sc1a2oup2", height = pList[input$sc1a2psz]) ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1a2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
    }) 
  output$sc1a2oup2.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
    }) 
  
  
  ### Plots for tab a3 
  output$sc1a3sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a3sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a3sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a3oup1 <- renderPlot({ 
    req(input$goButton_umap)
    scDRgene(sc1conf, sc1meta, isolate(input$sc1a3drX), isolate(input$sc1a3drY), isolate(input$sc1a3inp1),  
             isolate(input$sc1a3sub1), isolate(input$sc1a3sub2), 
             "sc1gexpr.h5", sc1gene, 
             isolate(input$sc1a3siz), isolate(input$sc1a3col1), isolate(input$sc1a3ord1), 
             isolate(input$sc1a3fsz), isolate(input$sc1a3asp), isolate(input$sc1a3txt)) 
  }) 
  output$sc1a3oup1.ui <- renderUI({ 
    withLoader(plotOutput("sc1a3oup1", height = pList[input$sc1a3psz]) ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1a3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   paste(input$sc1a3inp1,collapse = "_"),".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, isolate(input$sc1a3drX), isolate(input$sc1a3drY), isolate(input$sc1a3inp1),  
                      isolate(input$sc1a3sub1), isolate(input$sc1a3sub2), 
                      "sc1gexpr.h5", sc1gene, 
                      isolate(input$sc1a3siz), isolate(input$sc1a3col1), isolate(input$sc1a3ord1), 
                      isolate(input$sc1a3fsz), isolate(input$sc1a3asp), isolate(input$sc1a3txt)) ) 
    }) 
  output$sc1a3oup1.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   paste(input$sc1a3inp1,collapse = "_"),".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, 
      plot = scDRgene(sc1conf, sc1meta, isolate(input$sc1a3drX), isolate(input$sc1a3drY), isolate(input$sc1a3inp1),  
                      isolate(input$sc1a3sub1), isolate(input$sc1a3sub2), 
                      "sc1gexpr.h5", sc1gene, 
                      isolate(input$sc1a3siz), isolate(input$sc1a3col1), isolate(input$sc1a3ord1), 
                      isolate(input$sc1a3fsz), isolate(input$sc1a3asp), isolate(input$sc1a3txt)) ) 
    }) 

  ### Plots for tab b2 
  output$sc1b2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1b2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1b2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1b2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1b2oup1 <- renderPlot({ 
    scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,   
             input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
             input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) 
  }) 
  output$sc1b2oup1.ui <- renderUI({ 
    withLoader(plotOutput("sc1b2oup1", height = pList2[input$sc1b2psz]) ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1b2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
    }) 
  output$sc1b2oup1.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
    }) 
  output$sc1b2oup2 <- renderPlot({ 
    scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) 
  }) 
  output$sc1b2oup2.ui <- renderUI({ 
    withLoader(plotOutput("sc1b2oup2", height = "300px") ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1b2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
    }) 
  output$sc1b2oup2.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
    }) 
  output$sc1b2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(sc1conf, sc1meta, input$sc1b2inp1, input$sc1b2inp2, 
                         input$sc1b2sub1, input$sc1b2sub2, "sc1gexpr.h5", sc1gene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
  
  
  ### Plots for tab c1 
  output$sc1c1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  
  output$sc1c1oup <- renderPlot({
    req(input$goButton_vio)
    scVioBox(sc1conf, sc1meta, isolate(input$sc1c1inp1), isolate(input$sc1c1inp2),
             isolate(input$sc1c1sub1), isolate(input$sc1c1sub2),
             "sc1gexpr.h5", sc1gene, isolate(input$sc1c1typ), isolate(input$sc1c1pts),
             isolate(input$sc1c1siz), isolate(input$sc1c1fsz))
  })
  output$sc1c1oup.ui <- renderUI({ 
    withLoader(plotOutput("sc1c1oup", height = pList2[input$sc1c1psz]) ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1c1oup.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   paste(input$sc1c1inp2,collapse = "_"),".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1c1oup.h, width = input$sc1c1oup.w, useDingbats = FALSE, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, isolate(input$sc1c1inp2), 
                      input$sc1c1sub1, input$sc1c1sub2, 
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
    }) 
  output$sc1c1oup.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   paste(input$sc1c1inp2,collapse = "_"),".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1c1oup.h, width = input$sc1c1oup.w, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, isolate(input$sc1c1inp2), 
                      input$sc1c1sub1, input$sc1c1sub2, 
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
    }) 
  
  
  ### Plots for tab c2 
  output$sc1c2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1c2oup <- renderPlot({ 
    scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
           input$sc1c2sub1, input$sc1c2sub2, 
           input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) 
  }) 
  output$sc1c2oup.ui <- renderUI({ 
    withLoader(plotOutput("sc1c2oup", height = pList2[input$sc1c2psz]) ,type = "html",loader = "dnaspin") 
  }) 
  output$sc1c2oup.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                   input$sc1c2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1c2oup.h, width = input$sc1c2oup.w, useDingbats = FALSE, 
      plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                    input$sc1c2sub1, input$sc1c2sub2, 
                    input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
    }) 
  output$sc1c2oup.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                   input$sc1c2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1c2oup.h, width = input$sc1c2oup.w, 
      plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                    input$sc1c2sub1, input$sc1c2sub2, 
                    input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
    })
  
  ######## Plots for tab of gene signature tab expression ######
  output$sc1d1_list.ui <- renderUI({
    selected_set <- input$sc1d1_set 
    textAreaInput("sc1d1inp", HTML("Insert/Remove genes from list"), 
                  height = "400px", 
                  value = paste0(msigdbr_list[[selected_set]],collapse = "\n")) %>% 
      helper(type = "inline", size = "m", fade = TRUE, 
             title = "List of genes to plot on gene signature", 
             content = c("Input genes to plot", 
                         "- Genes should be separated by comma, semicolon or newline"))
  }) 
  

  output$sc1d1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1d1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1d1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1d1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1d1oupTxt <- renderUI({ 
    geneList = scGeneList(input$sc1d1inp, sc1gene) 
    if(nrow(geneList) == 0){ 
      HTML("No gene list was provided!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes found in the dataset") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  
  #Change scBubbHeat to function of interest
  output$sc1d1oup <- renderPlot({
    req(input$goButton_exp)
    #Only plot after using goButton, thus use isolate
    plot_gene_signature(sc1conf, sc1meta, isolate(input$sc1d1inp), isolate(input$sc1d1grp), isolate(input$sc1d1plt),
                        isolate(input$sc1d1sub1), isolate(input$sc1d1sub2), "sc1gexpr.h5", sc1gene,
                        isolate(input$sc1d1cols), isolate(input$sc1d1fsz),isolate(input$sc1d1col1),
                        isolate(input$sc1d1siz),isolate(input$sc1d1lab1),"expression")
  })
  
  output$sc1d1oup.ui <- renderUI({ 
    withLoader(plotOutput("sc1d1oup", height = pList3[input$sc1d1psz]),type = "html",loader = "dnaspin") 
  }) 
  output$sc1d1oup.pdf <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1d1plt,"_",input$sc1d1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = last_generated_bar ) 
    }) 
  output$sc1d1oup.png <- downloadHandler( 
    filename = function() { paste0("Tcells_",input$sc1d1plt,"_",input$sc1d1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = last_generated_bar )
    })
  
  ######## Plots for tab of gene signature tab correlation ######
  output$sc1d2_list.ui <- renderUI({
    selected_set <- input$sc1d2_set 
    textAreaInput("sc1d2inp", HTML("Insert/Remove genes from list"), 
                  height = "400px", 
                  value = paste0(msigdbr_list[[selected_set]],collapse = "\n")) %>% 
    helper(type = "inline", size = "m", fade = TRUE, 
           title = "List of genes to plot on gene signature", 
           content = c("Input genes to plot", 
                       "- Genes should be separated by comma, semicolon or newline"))
  }) 
  
  output$sc1d2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1d2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1d2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1d2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1d2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1d2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1d2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1d2oupTxt <- renderUI({ 
    geneList = scGeneList(input$sc1d2inp, sc1gene) 
    if(nrow(geneList) == 0){ 
      HTML("No gene list was provided!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes found in the dataset") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  
  #Change scBubbHeat to function of interest
  output$sc1d2oup <- renderPlot({
    req(input$goButton_correlation)
    #Only plot after using goButton, thus use isolate
    plot_gene_signature(sc1conf, sc1meta, isolate(input$sc1d2inp), isolate(input$sc1d2grp), isolate(input$sc1d2plt),
                        isolate(input$sc1d2sub1), isolate(input$sc1d2sub2), "sc1gexpr.h5", sc1gene,
                        isolate(input$sc1d2cols), isolate(input$sc1d2fsz),isolate(input$sc1d2col1),
                        isolate(input$sc1d2siz),isolate(input$sc1d2lab1),"correlation")
  })
  
  output$sc1d2oup.ui <- renderUI({ 
    withLoader(plotOutput("sc1d2oup", height = pList3[input$sc1d2psz]),type = "html",loader = "dnaspin") 
  }) 
  output$sc1d2oup.pdf <- downloadHandler( 
    filename = function() { paste0("Correlation_plot.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1d2oup.h, width = input$sc1d2oup.w, 
      plot = last_generated ) 
    }) 
  output$sc1d2oup.png <- downloadHandler( 
    filename = function() { paste0("Correlation_plot.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1d2oup.h, width = input$sc1d2oup.w, 
      plot = last_generated )
    })
  
}) 









