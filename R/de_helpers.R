#' Create edgeR object
#'
#' @title
#' @import edgeR
#' @return
#' @author dylanmr
#' @export
edger_prep <- function(seur, trt_group = "geno", lib.size = 1e4,
                       filter_column = NULL, filter_value = NULL,
                       assay_name = "RNA", feature_pattern_exclude = "^Gm|Rik$", 
                       dataset = NULL, celltypes = celltypes, celltype_column = NULL,
                       min.count = 1, min.total.count = 5) {
  
  # Validate inputs
  stopifnot("Seurat" %in% class(seur), assay_name %in% names(seur@assays))
  DefaultAssay(seur) <- "RNA"
  
  # Subset cells based on metadata column and filter features
  cells_to_keep <- colnames(seur)[grepl(paste0("^",celltypes, "$"), seur[[celltype_column]][,1])]
  print(length(cells_to_keep))
  features_to_keep <- rownames(seur@assays[[assay_name]])[!grepl(feature_pattern_exclude, rownames(seur@assays[[assay_name]]))]
  seur <- subset(seur, cells = cells_to_keep, features = features_to_keep)
  
  if(!is.null(dataset)){
    cells_to_keep <- colnames(seur)[grepl(paste0("^",dataset, "$"), seur[["dataset"]][,1])]
    print(length(cells_to_keep))
    seur <- subset(seur, cells = cells_to_keep)
  }
  
  if(!is.null(filter_column)){
    cells_to_keep <- colnames(seur)[grepl(paste0("^",filter_value, "$"), seur[[filter_column]][,1])]
    print(length(cells_to_keep))
    seur <- subset(seur, cells = cells_to_keep)
  }
  
  y <- Seurat2PB(seur, sample = "hash.mcl.ID", cluster= trt_group)
  keep.samples <- y$samples$lib.size > lib.size
  y <- y[, keep.samples]
  meta <- seur[[]] %>% distinct(hash.mcl.ID, .keep_all = T) 
  meta <- meta[match(y$samples$sample, meta$hash.mcl.ID),]
  y$samples <- bind_cols(meta, y$samples)
  
  keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count = min.count, min.total.count = min.total.count)
  table(keep.genes)
  y <- y[keep.genes, , keep=FALSE]
  y <- normLibSizes(y)
  
  return(y)
  
}
