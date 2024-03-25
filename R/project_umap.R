#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
project_umap <- function(ref, query, label_to_transfer="predicted.id") {
  
  anchors <- FindTransferAnchors(reference = ref, query = query, reference.assay = "integrated", query.assay = "integrated", 
                      dims = 1:30, reference.reduction = "pca")

  query <- MapQuery(anchorset = anchors, reference = ref, query = query,
                             refdata = list(celltype = label_to_transfer), reference.reduction = "pca", reduction.model = "umap")
  
  return(query)

}
