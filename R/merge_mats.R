#' .. content for \description{} (no empty lines) ..
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
