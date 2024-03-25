#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
gen_xenium <- function(path, min_umi = 3) {

  obj <- LoadXenium(path, fov="fov")
  pattern <- "__(R.*?)__"
  obj[["section"]] <- str_match_all(path, pattern)[[1]][,2]
  pattern <- "__(00.*?)__"
  obj[["slide"]] <- str_match_all(path, pattern)[[1]][,2]
  obj <- subset(obj, subset = nCount_Xenium > min_umi)
  return(obj)

}
