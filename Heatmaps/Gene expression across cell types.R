#### Heatmap of normalized marker gene expression across shared major cell types ####

# Codes created by Sam Vanmassenhove

# As an example, we present the expression of marker genes across major cell types.
# The same approach has also been applied at the cell subtype level.

### Load the library
library(Seurat) 

### Run the function 
cleanHeatmap3 <- function(object, 
                          genes = NULL,
                          assay = "RNA",
                          group.by = "SampleID",
                          idents = NULL,
                          invert = F,
                          space.width = 0.92,
                          space.height = 0.92,
                          facet = TRUE,
                          xlab=NULL,
                          ylab=NULL,
                          line.draw=T,
                          xticks=T,
                          yticks=T,
                          tick.angle=90,
                          axis.text.size=14,
                          facet.text.size=14, 
                          slot="data") {
  # Get the standard list of genes if none provided
  if (is.null(genes)) {
    genes <- list(
      "Tcell" = c("CD3D", "CD3E", "CD8A", "CD4"),
      "NK" = c("NCAM1"),
      "Bcell" = c("CD79A"),
      "Myeloid" = c("CD68", "LYZ", "AIF1"),
      "Mast" = c("MS4A2", "TPSAB1", "CPA3"),
      "pDC" = c("LILRA4", "CXCR3"),
      "cDC" = c("CLEC9A", "XCR1", "CD1C", "CCR7", "CCL17", "CCL19"),
      "EC" = c("CLDN5", "PECAM1", "VWF"),
      "Fibroblast" = c("COL1A1", "COL10A1", "COL4A1", "BGN", "DCN"))
  }
  DefaultAssay(object) <- assay
  Idents(object) <- object@meta.data[[group.by]]
  markers <- as.vector(unlist(genes))
  markers.df <- data.frame(Gene = markers)
  # Handle metadata non-gene "markers" (such as nCount_RNA)
  meta.df <- NULL
  meta.markers <- markers[markers %in% colnames(object@meta.data)]
  if (length(meta.markers) != 0) {
    meta.df <- aggregate(object[[meta.markers]], list(Idents(object)), mean)
    meta.df <- data.table::transpose(meta.df, make.names = "Group.1")
    rownames(meta.df) <- meta.markers
  }
  # Handle gene markers
  markers <- markers[markers %in% rownames(object)]
  expr.df <- AverageExpression(object, 
                               assays = assay, 
                               features = markers, 
                               group.by = group.by,
                               slot = slot)[[assay]]
  expr.df <- as.data.frame(expr.df)
  # Join the two types
  all.markers <- c(markers, meta.markers)
  expr.df <- rbind(expr.df, meta.df)
  # Join the markers with the expression dataframe
  expr.df <- expr.df %>%
    tibble::rownames_to_column(var = "Gene")
  # Continue with the requested idents only
   if (is.null(idents) & !invert) idents <- as.vector(unique(Idents(object)))
  expr.df <- if (!invert) dplyr::select(expr.df, c("Gene", idents)) else dplyr::select(expr.df, -idents)
  # Calculate the z-score
  expr.df <- expr.df %>%
    tidyr::gather(key = cell.type , value = expr, -Gene) %>%
    group_by(Gene) %>%
    dplyr::mutate(zscore = (expr - mean(expr))/(sd(expr) + 0.01))
  # Add subtype column to create facets if requested
  if (facet & class(genes) == "list") {
    expr.df$subtype <- NULL
    for (subtype in names(genes)) {
      expr.df[expr.df$Gene %in% genes[[subtype]], "subtype"] = subtype
    }
    # Order facets according to listed order
    expr.df$subtype <- factor(expr.df$subtype, ordered=TRUE, levels=names(genes))
  }
  # Order the genes according to listed order
  expr.df$Gene <- factor(expr.df$Gene, ordered=TRUE, levels=rev(unique(all.markers)))
  # Reorder the x-axis idents if they can be cast to numeric: Seurat sometimes orders clusters alphabetically rather than numerically
  if (suppressWarnings(!any(is.na(as.numeric(idents))))) {
    ordered.levels = sort(unique(as.numeric(idents)))
    expr.df$cell.type <- factor(expr.df$cell.type, ordered=T, levels=ordered.levels)
  } 
  # Reorder the cell type idents to their original order if it was already factorized
  if (is.factor(Idents(object))) { 
    expr.df$cell.type <- factor(expr.df$cell.type, levels = if (!is.factor(object@meta.data[[group.by]])) unique(expr.df$cell.type) else levels(object@meta.data[[group.by]]))
    # Drop empty levels if they exist
    expr.df$cell.type <- droplevels(expr.df$cell.type)
  }
  # Create the heatmap =========================================================
  heatmap <- ggplot(data = expr.df) +
    geom_tile(aes(x = cell.type, y = Gene, fill = zscore, width=space.width, height=space.height)) +
    xlab(label = xlab) +
    ylab(label = ylab) + 
    scale_fill_gradient2(low = "darkblue", mid= "white", high = "darkred")
  # Add facets if requested
  if (facet & class(genes) == "list") {
    heatmap <- heatmap + facet_grid(rows=vars(subtype), scales="free_y", space="free_y")
  }
  # Make horizontal facet labels with a publishable theming
  p <- heatmap + theme(strip.text.y = element_text(angle=0),
                       axis.text.x = element_text(angle=tick.angle, hjust=1, size=axis.text.size),
                       axis.text.y = element_text(size=axis.text.size),
                       strip.text = element_text(size=facet.text.size),
                       panel.background = element_rect(fill=NA))
  # Draw the panel lines if required
  if (line.draw) {
    p <- p + theme(strip.background = element_rect(fill=NA),
                   panel.background = element_rect(fill=NA, color="black", linetype="solid"),
                   panel.border = element_rect(fill=NA))
  }
  # Remove the axis ticks if required
  if (!xticks) p <- p + theme(axis.ticks.x = element_line(colour = NA))
  if (!yticks) p <- p + theme(axis.ticks.y = element_line(colour = NA))
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
heatmap <- cleanHeatmap3(object, list("T cell"=c("CD3D", "CD3E", "TRAC"), "NK"= c("NCAM1", "KLRC3", "KLRC1"),"B cell"= c("CD79A", "CD79B", "MS4A1"),"Macro/Mono" = c("CD68", "LYZ", "AIF1"),
                                  "DC" = c("LILRA4", "XCR1", "CLEC9A"), "Mast cell"=c("MS4A2", "TPSAB1", "CPA3"), "Cancer/epit cell" = c("EPCAM", "KRT7", "KRT18"), "Fibroblast"=c("COL1A1", "COL10A1", "BGN"), 
                                  "EC"= c("CLDN5", "PECAM1", "VWF")), group.by = "Intermediatecelltype_annotation")
### Visualize the heatmap 
heatmap

### Save the heatmap
ggsave("/path/Heatmap_Intermediatecelltype_annotation.pdf", plot = heatmap, width = 5.5, height = 10) 



#### Analyze T-cell Reactivity score across T-cell subtypes ####
# Codes created by Sam Vanmassenhove

### Load libraries 
library(pheatmap)
library(Seurat)

### Read T/NK-cell object from individual cancer types (see file "Tcell analysis.R") 
Tcell <- readRDS("/path/Tcell_object.rds")

### Function to add expression data for multiple genes to metadata
add_gene_expression <- function(seurat_obj, gene_names) {
  for (gene_name in gene_names) {
    # Fetch expression data for the current gene for all cells
    gene_expression <- FetchData(seurat_obj, vars = gene_name)
    
    # Add the expression data as a new column in the metadata
    seurat_obj <- AddMetaData(seurat_obj, metadata = gene_expression, col.name = paste0(gene_name, "_expression"))
  }
  return(seurat_obj)
}

### List of genes for which you want to add expression data
gene_list <- c("ITGAE", "PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1", "TNFRSF9", "TNFRSF18")

### Call the function to add expression data for multiple genes
Tcell <- add_gene_expression(Tcell, gene_list)

### Print the updated metadata
head(Tcell@meta.data)

colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

### add annotations 
features = c('ITGAE_expression', 'PDCD1_expression', 'CTLA4_expression', 'LAG3_expression', 'GZMB_expression', 'PRF1_expression', 'TNFRSF9_expression', 'TNFRSF18_expression')

tf <- Tcell@meta.data[,c("Minorcelltype_annotation",features)] %>%
  group_by(Minorcelltype_annotation) %>% 
  summarise_all(~mean(.x, na.rm = TRUE)) %>%
  as.data.frame()

rownames(tf) <- tf$Minorcelltype_annotation
tf$Minorcelltype_annotation <- NULL

### Create heatmap
heatmap <- pheatmap(
  tf,
  color = colors,
  treeheight_col = 0,
  cluster_cols = FALSE,
  border_color = TRUE,
  cellwidth = 15,
  cellheight = 15,
  fontsize = 12,
  main = "T-cell reactivity gene score across T-cell subtypes",
  # annotation = annotation, 
  # annotation_colors = list(Signatures = Signatures),
  scale = "column", #scaling for columns
  angle_col = 90, # Rotate x-axis labels by 90 degrees
  labels_col = features  # Use modified features vector here
)

### Visualize the heatmap 
heatmap

### Save the heatmap
ggsave("/path/Heatmap_TcellsReactivity_genes_acrossTcellsubtypes.pdf", plot = heatmap) 
