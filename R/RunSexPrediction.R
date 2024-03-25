#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param x
#' @return
#' @author dylanmr
#' @export
RunSexPrediction <- function(x, features = c("Xist","Tsix"), assay= "SCT") {

  sex_genes <- AverageExpression(x, features = features, assays = assay, group.by = "hash.mcl.ID")$RNA
  x@meta.data <- full_join(x@meta.data, enframe(factor(Mclust(colMeans(sex_genes),G=2)$classification), 
                                                name = "hash.mcl.ID", value="sex_predicted"), by=c("hash.mcl.ID"))
  rownames(x@meta.data) <- colnames(x)
  return(x)

}
