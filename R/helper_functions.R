#' Convert seurat to h5 object
#' @title
#' @import zellkonverter
#' @return 
#' @author dylanmr
#' @export
convert_to_h5 <- function(obj, filename) {
  
  DefaultAssay(obj) <- "RNA"
  sce <- as.SingleCellExperiment(obj)
  colData(sce) <- colData(sce)[colSums(is.na(colData(sce))) == 0]
  writeH5AD(sce, file = filename)
  
}

#' Convert seurat to loom file
#' @return
#' @import SeuratDisk
#' @author dylanmr
#' @export
create_loom <- function(obj, filepath, subset_value, subset_column) {
  
  if(!is.null(subset_value)) {
    obj <- filter_cells_by_column(obj, subset_column, subset_value)
  }
  
  if(!is.null(obj[["integrated"]])) {
    obj[["integrated"]] <- NULL
  }
  
  DefaultAssay(obj) <- "RNA"
  loom_ma <- SeuratDisk::as.loom(obj, filename = filepath, verbose = FALSE)
  loom_ma$close_all()
  
}

#' code to merge matrices
#'
#' .. content for \details{} ..
#'
#' @title
#' @param filtered_seurat
#' @return
#' @author dylanmr
#' @export
merge_mats <- function(filtered_seurat) {
  
  mat_list <- purrr::map(filtered_seurat, function(x) {
    cts <- x@assays$RNA@counts
    colnames(cts) <- paste0(colnames(cts), "_", x$orig.ident)
    return(cts)
  })  
  
  meta_df <- purrr::map_dfr(filtered_seurat, function(x) {
    meta <- x@meta.data
    rownames(meta) <- paste0(rownames(meta), "_", x$orig.ident)
    return(meta)
  })
  
  merged_mat <- RowMergeSparseMatrices(mat_list[[1]], mat_list[-1])
  seur <- CreateSeuratObject(merged_mat, meta.data = meta_df)
  
  return(seur)
  
}

#' Helper to create new column from existing columns in seurat object
#'
#' @title
#' @return
#' @author dylanmr
#' @export
create_new_column_seurat <- function(obj, selected_columns = c("treatment", "time", "geno"), new_column_name = "group") {
  
  df <- obj[[]]
  
  # Check if all selected columns are in the dataframe
  if(!all(selected_columns %in% names(df))) {
    stop("Not all selected columns are in the dataframe")
  }
  
  # Create the new column by pasting the selected columns together
  # The 'apply' function is used to iterate over rows
  obj[[new_column_name]] <- apply(df[, selected_columns, drop = FALSE], 1, paste, collapse = "_")
  
  # Return the modified dataframe
  return(obj)
  
}
