#### Heatmap of normalized marker gene expression across shared major cell types per cancer type ####

# Codes created by Sam Vanmassenhove

# As an example, we present the expression of marker genes across major cell types.
# The same approach has also been applied at the cell subtype level.

# Load the library
library(Seurat) 

# Run the function 
Heatmap_tumortypes <- function (object, 
                                genes = NULL, 
                                assay = "RNA", 
                                group.by = "SampleID", 
                                order.tumortype = NULL,
                                order.celltype = NULL,
                                idents = NULL, 
                                invert = F, 
                                space.width = 0.92, 
                                space.height = 0.92, 
                                facet = TRUE, 
                                xlab = NULL, 
                                ylab = NULL, 
                                line.draw = T, 
                                xticks = T, 
                                yticks = T, 
                                tick.angle = 90, 
                                axis.text.size = 14, 
                                facet.text.size = 14, 
                                slot = "data") 
{
  
  DefaultAssay(object) <- assay
  markers <- as.vector(unlist(genes))
  #markers.df <- data.frame(Gene = markers)
  meta.df <- NULL
  markers <- markers[markers %in% rownames(object)]
  
  # Clean names to remove forbidden character "_"
  orig.names1 = unique(object@meta.data[[group.by[1]]]) %>% as.character()
  orig.names2 = unique(object@meta.data[[group.by[2]]]) %>% as.character()
  object@meta.data[[group.by[1]]] = as.character(object@meta.data[[group.by[1]]])
  object@meta.data[[group.by[2]]] = as.character(object@meta.data[[group.by[2]]])
  object@meta.data[[group.by[1]]] = factor(object@meta.data[[group.by[1]]], 
                                           levels = orig.names1,
                                           labels = gsub(pattern="_", x=orig.names1, replacement = ";"))
  object@meta.data[[group.by[2]]] = factor(object@meta.data[[group.by[2]]], 
                                           levels = orig.names2,
                                           labels = gsub(pattern="_", x=orig.names2, replacement = ";"))
  object@meta.data[[group.by[1]]] = as.character(object@meta.data[[group.by[1]]])
  object@meta.data[[group.by[2]]] = as.character(object@meta.data[[group.by[2]]])
  
  
  # Use AverageExpression function to calculate the mean per group
  expr.df <- AverageExpression(object, assays = assay, features = markers, 
                               group.by = group.by, slot = slot)[[assay]] %>% as.data.frame %>% rownames_to_column("Gene")
  expr.df <- as.data.frame(expr.df) %>%
    pivot_longer(names_to = group.by, values_to = "score", names_sep = "_", cols = -Gene)
  expr.df <- expr.df %>% group_by(Gene, !!sym(group.by[2])) %>% 
    dplyr::mutate(zscore = (score - mean(score))/(sd(score) + 0.01))
  
  # Fix the names again
  expr.df[[group.by[1]]] = gsub(pattern = ";", expr.df[[group.by[1]]], replacement="_")
  expr.df[[group.by[2]]] = gsub(pattern = ";", expr.df[[group.by[2]]], replacement="_")
  
  # Set the levels for the order of celltypes, tumor types and genes
  expr.df$Gene = factor(expr.df$Gene, levels = genes)
  expr.df[[group.by[1]]] = factor(expr.df[[group.by[1]]], levels = if (is.null(order.celltype)) orig.names1 else order.celltype)
  expr.df[[group.by[2]]] = factor(expr.df[[group.by[2]]], levels = if (is.null(order.tumortype)) orig.names2 else order.tumortype)
  
  heatmap <- ggplot(data = expr.df) + 
    geom_tile(aes(x = !!sym(group.by[1]), y = !!sym(group.by[2]), fill = zscore, 
                  width = space.width, height = space.height)) + 
    xlab(label = xlab) + ylab(label = ylab) + 
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") + 
    facet_grid(rows = "Gene", scales = "free", space = "free")
  
  
  p <- heatmap + 
    theme(strip.text.y = element_text(angle = 0), 
          axis.text.x = element_text(angle = tick.angle, hjust = 1, vjust = 0.5, color = "black", 
                                     size = axis.text.size), 
          axis.text.y = element_text(size = axis.text.size, color = "black"), 
          strip.text = element_text(size = facet.text.size), panel.background = element_rect(fill = NA))
  
  if (line.draw) {
    p <- p + theme(strip.background = element_rect(fill = NA), 
                   panel.background = element_rect(fill = NA, color = "black", 
                                                   linetype = "solid"), panel.border = element_rect(fill = NA))
  }
  
  if (!xticks) 
    p <- p + theme(axis.ticks.x = element_line(colour = NA))
  if (!yticks) 
    p <- p + theme(axis.ticks.y = element_line(colour = NA))
  return(p)
}



###### Accessing high-quality single-cell data
# To generate a Seurat object containing all shared pan-cancer high-quality cells, users can access the read count data for each individual cancer type on our lab’s website (https://lambrechtslab.sites.vib.be/en/dataaccess).
# After performing annotation analysis (refer to "Individual Cancer Types Analysis" and the cell subtype annotations in the "Cell Subtype Annotation" folder), users can explore the metadata of our pan-cancer atlas, available as "Lodietall_metadata.csv" within the "Master_files" folder. 
# This metadata file includes both general sample information and detailed cell subtype annotations for further reference.
# Once this object has been generated, proceed with the following analysis. 

### Read object with all cells (no doublets and no low quality cells)
object <- readRDS("/path/object.RDS")

### Create the heatmap
heatmap <- Heatmap_tumortypes(object, genes = c("PDCD1", "CD274", "PDCD1LG2"), group.by=c("Intermediatecelltype_annotation", "TumorType"), 
                          order.tumortype = rev(c("BC_early", "BC_adv", "CC", "CRC", "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv")), 
                          order.celltype = c("T cell", "NK",  "B cell", "Macro/Mono", "DC", "Mast cell", "Cancer/epit cell","Fibroblast","EC"))

# Visualize the heatmap 
heatmap

# Save the heatmap
ggsave("/path/Heatmap_Intermediatecelltype_acrossCancerTypes.pdf", plot = heatmap, width = 5, height = 8.5) 