print.noquote("bboeckx, check reults carefully")
print.noquote("load also Library (Matrix)") 

ReadEnsemble <- function (data.dir = NULL, gene.column = 1, unique.features = TRUE, 
    strip.suffix = FALSE) 
{
    full.data <- list()
    for (i in seq_along(along.with = data.dir)) {
        run <- data.dir[i]
        if (!dir.exists(paths = run)) {
            stop("Directory provided does not exist")
        }
        barcode.loc <- file.path(run, "barcodes.tsv")
        gene.loc <- file.path(run, "genes.tsv")
        features.loc <- file.path(run, "features.tsv.gz")
        matrix.loc <- file.path(run, "matrix.mtx")
        pre_ver_3 <- file.exists(gene.loc)
        if (!pre_ver_3) {
            addgz <- function(s) {
                return(paste0(s, ".gz"))
            }
            barcode.loc <- addgz(s = barcode.loc)
            matrix.loc <- addgz(s = matrix.loc)
        }
        if (!file.exists(barcode.loc)) {
            stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
        }
        if (!pre_ver_3 && !file.exists(features.loc)) {
            stop("Gene name or features file missing. Expecting ", 
                basename(path = features.loc))
        }
        if (!file.exists(matrix.loc)) {
            stop("Expression matrix file missing. Expecting ", 
                basename(path = matrix.loc))
        }
        data <- readMM(file = matrix.loc)
        cell.names <- readLines(barcode.loc)
        if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
            cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                FUN = ExtractField, field = 1, delim = "-")))
        }
        if (is.null(x = names(x = data.dir))) {
            if (i < 2) {
                colnames(x = data) <- cell.names
            }
            else {
                colnames(x = data) <- paste0(i, "_", cell.names)
            }
        }
        else {
            colnames(x = data) <- paste0(names(x = data.dir)[i], 
                "_", cell.names)
        }
        feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
            yes = gene.loc, no = features.loc), header = FALSE, 
            stringsAsFactors = FALSE)
        if (any(is.na(x = feature.names[, gene.column]))) {
            warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
                call. = FALSE, immediate. = TRUE)
            na.features <- which(x = is.na(x = feature.names[, 
                gene.column]))
            replacement.column <- ifelse(test = gene.column == 
                2, yes = 1, no = 2)
            feature.names[na.features, gene.column] <- feature.names[na.features, 
                replacement.column]
        }
        if (unique.features) {
            fcols = ncol(x = feature.names)
            if (fcols < gene.column) {
                stop(paste0("gene.column was set to ", gene.column, 
                  " but feature.tsv.gz (or genes.tsv) only has ", 
                  fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                  fcols, "."))
            }
            rownames(x = data) <- make.unique(names = feature.names[, 
                gene.column])
        }
        if (ncol(x = feature.names) > 2) {
            data_types <- factor(x = feature.names$V3)
            lvls <- levels(x = data_types)
            if (length(x = lvls) > 1 && length(x = full.data) == 
                0) {
                message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
            }
            expr_name <- "Gene Expression"
            if (expr_name %in% lvls) {
                lvls <- c(expr_name, lvls[-which(x = lvls == 
                  expr_name)])
            }
            data <- lapply(X = lvls, FUN = function(l) {
                return(data[data_types == l, , drop = FALSE])
            })
            names(x = data) <- lvls
        }
        else {
            data <- list(data)
        }
        full.data[[length(x = full.data) + 1]] <- data
    }
    list_of_data <- list()
    for (j in 1:length(x = full.data[[1]])) {
        list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
            FUN = `[[`, j))
        list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
    }
    names(x = list_of_data) <- names(x = full.data[[1]])
    if (length(x = list_of_data) == 1) {
        return(list_of_data[[1]])
    }
    else {
        return(list_of_data)
    }
}
