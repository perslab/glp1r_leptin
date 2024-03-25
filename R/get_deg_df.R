#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
get_deg_df <- function(dds, contrast = NULL) {

  if(!is.null(contrast)) {
    res <- 
      results(dds) %>% 
      data.frame() %>% 
      rownames_to_column("gene")
    res <- 
      res %>% 
      mutate(celltype = rep(unique(dds$labels), dim(res)[1]))
  }
  else {
    res <- 
      results(dds) %>% 
      data.frame() %>% 
      rownames_to_column("gene")
    res <- 
      res %>% 
      mutate(celltype = rep(unique(dds$final_labels), dim(res)[1]))
  }
  
  return(res)

}
