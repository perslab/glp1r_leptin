#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
id_enriched_clusters <- function(obj, cluster_column, sample_column, grouping_factor) {
  
  meta <- obj@meta.data
  meta$group <- ifelse(grepl("Sun|Fl", meta$dataset, ignore.case = T), "mm","dr") %>%  as.factor
  res <- speckle::propeller(clusters = meta[[cluster_column]], sample = meta[[sample_column]], group = meta[[grouping_factor]])
  ids <- res %>% dplyr::filter(sign(Tstatistic)<0) %>%  filter(FDR<0.01) %>% mutate(diff = log2(PropMean.dr/PropMean.mm)) %>%  filter(diff< (-1)) %>%  pull(1)
  return(ids)
  
}
