#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param obj
#' @param out
#' @return
#' @author dylanmr
#' @export
create_loom <- function(obj, filepath) {

  DefaultAssay(obj) <- "RNA"
  loom_ma <- as.loom(obj, filename = filepath, verbose = FALSE)
  loom_ma$close_all()

}
