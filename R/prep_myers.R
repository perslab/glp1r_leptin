#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
prep_myers <- function(path, min.features, min.cells, name = "leprgfp") {

  obj <- Read10X(here::here(path), gene.column = 1, strip.suffix = T)
  colnames(obj) <- stri_sub(stri_sub(colnames(obj),to =  -2), 2)
  meta <- read.csv(here::here(path, "metadata.csv.gz"), row.names = "ID")
  obj <- CreateSeuratObject(obj, min.features = min.features, min.cells = min.cells, project = name, meta.data = meta)
  obj <- subset(obj, subset = Cluster_region == "HY")
  obj <- subset(obj, subset = Lepr_cluster %in% c("L10.GABA","L8.Nts"), invert=T)
  clustokeep <- obj$Lepr_cluster[grep("^L", obj$Lepr_cluster)] %>% unique()
  obj <- subset(obj, subset = Lepr_cluster %in% clustokeep)
  obj <- process_seurat(obj, method="log", res=0.8, dims=30)
  return(obj)

}
