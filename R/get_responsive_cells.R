#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param obwt_cacoa
#' @return
#' @author dylanmr
#' @export
get_responsive_cells <- function(x, cutoff = 0.05, truth = NULL) {

  pvals <- x$test.results$expression.shifts$dists.per.type %>% 
    enframe() %>% 
    mutate(pv = x$test.results$expression.shifts$pvalues,
           padj = x$test.results$expression.shifts$padjust) %>% 
    unnest(value)
  
  impacted_pops <- 
    pvals %>%  
    filter(padj < cutoff) %>% 
    pull(name) %>% 
    unique()
  
  impacted_pops <- ifelse(grepl(".", impacted_pops), gsub("[.].*", "", impacted_pops) %>%  unique(), unique(impacted_pops))

  return(impacted_pops)

}
