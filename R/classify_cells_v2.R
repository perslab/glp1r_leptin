#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export

classify_cells_v2 <- function(obj, path =  here::here("cell_markers/coarse_markers.txt"), 
                           clustering_res = NULL,
                           clustering_dims = NULL) {
  
  cm <- obj@assays$RNA@counts
  cm_norm <- Matrix::t(obj@assays$RNA@data)
  emb <- obj@reductions$umap@cell.embeddings
  clusters <- setNames(obj@meta.data$seurat_clusters, rownames(obj@meta.data))
  clf_data <- getClassificationData(cm, markers = path)
  ann_by_level <- assignCellsByScores(graph = NULL, clf_data, clusters=clusters)
  obj[["ann_neuron_v2"]] <- ann_by_level$annotation$l1
  return(obj)
  
}
