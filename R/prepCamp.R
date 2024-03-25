#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param nameme1
#' @return
#' @author dylanmr
#' @export

prepCamp <- function(method) {
  
  camp <- scRNAseq::CampbellBrainData()
  seur <- process_seurat(camp, method = method, cluster=F, type="sce")
  seur <- AddMetaData(seur, metadata = data.frame(colData(camp)))
  return(seur)
  
}