#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export

classify_cells <- function(obj, path =  here::here("cell_markers/coarse_markers.txt"), 
                           clustering_res = clustering_res,
                           clustering_dims = clustering_dims) {
  
  obj <- process_seurat(obj, method = "log", res = clustering_res, dims = clustering_dims) 
  cm <- obj@assays$RNA@counts
  clusters <- setNames(obj@meta.data$seurat_clusters, rownames(obj@meta.data))
  clf_data <- getClassificationData(cm, markers = path)
  ann_by_level <- assignCellsByScores(graph = NULL, clf_data, clusters=clusters)
  obj[["ann_neuron"]] <- ann_by_level$annotation$l1
  return(obj)
  
}
