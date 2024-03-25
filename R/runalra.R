#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export

runalra <- function(obj) {

  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj)
  obj <- RunALRA(object = obj, use.mkl = F, slot = "data")

}
