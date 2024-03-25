#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param lep_responders
#' @param nameme1
#' @param nameme2
#' @param nameme3
#' @param nameme4
#' @return
#' @author dylanmr
#' @export

runspls <- function(obj, cells, assay, predictor, batch, scale=scale) {
  
  obj <- subset(obj, cells=cells)
  if(assay=="SCT") {
    obj <- process_seurat(obj, method = "glm", cluster=F, batch="hash_pool")
  }
  
  X <- GetAssayData(obj, slot = "data", assay=assay)
  X <- t(X[!grepl("^Gm|Rik$", rownames(X)),])
  Y <- factor(obj[[predictor]][,1]) 
  
  basic.plsda.model <- mint.splsda(X, Y, study = factor(obj[[batch]][,1]), ncomp = 2, keepX = 50, scale = scale) 
  gene_importance <- selectVar(basic.plsda.model)$value %>% arrange((value.var)) %>% rownames_to_column("gene")
  return(gene_importance)
}
