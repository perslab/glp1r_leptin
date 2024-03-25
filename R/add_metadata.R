#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param .
#' @param joinby
#' @return
#' @author dylanmr
#' @export

add_metadata <- function(obj, file = "sc_preprocessing/sequenced_mice_data_swapped_w_diet_ce.xlsx", joinby = "Hash") {
  
  data <- readxl::read_xlsx(file)  %>%  janitor::clean_names() 
  obj@meta.data <- left_join(obj@meta.data, data, by = c("hash.mcl.ID" = {{ joinby }}))
  rownames(obj@meta.data) <- colnames(obj)
  return(obj)
  
}