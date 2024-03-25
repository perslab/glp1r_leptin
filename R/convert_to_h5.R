#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
convert_to_h5 <- function(obj, filename) {

  DefaultAssay(obj) <- "RNA"
  sce <- as.SingleCellExperiment(obj)
  colData(sce) <- colData(sce)[colSums(is.na(colData(sce))) == 0]
  writeH5AD(sce, file = filename)

}
