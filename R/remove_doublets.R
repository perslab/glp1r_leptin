#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export

remove_doublets <- function(obj, k=30) {
  
  # Remove negative cells
  obj <- subset(obj, subset = HTO_mcl_classification.global == "Negative", invert=T)
  
  # Collect Data for scDblFinder
  dat <- GetAssayData(obj, slot = "data")
  dat <- dat[rownames(dat) %in% obj@assays$RNA@var.features,]
  doublets <- as.logical(ifelse(obj$HTO_mcl_classification.global == "Doublet", T, F))
  res <- scDblFinder::recoverDoublets(dat, doublets=doublets, samples=table(obj$hash.mcl.ID), k=k)
  print(paste(sum(res$predicted),"intra-sample doublets detected"))
  
  # Remove all doublets
  obj <- subset(obj, cells = Cells(obj)[which(res$predicted==T)], invert=T)
  obj <- subset(obj, subset= HTO_mcl_classification.global == "Singlet")
  
  return(obj)
  
}
