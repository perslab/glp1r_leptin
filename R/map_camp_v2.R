#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param obj
#' @param type
#' @param dims
#' @return
#' @author dylanmr
#' @export
map_camp_v2 <- function(obj = integrated_seurat, type = "coarse", dims = 30, hypomap) {

  if(type == "neuron") {
    camp <- readRDS("/projects/dylan/leptin_paper/data/Arc_Neurons_33_Clusters.rds")
    DefaultAssay(camp) <- "integrated"
    DefaultAssay(obj) <- "integrated"
  } else if(type == "coarse") {
    camp <- readRDS("/projects/dylan/leptin_paper/data/Arcuate_Full_Atlas.rds")
    DefaultAssay(camp) <- "RNA"
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj) %>% FindVariableFeatures()
  } else if(type == "hypomap") {
    hypomap <- hypomap
    obj <- NormalizeData(obj) %>% FindVariableFeatures()
  }
  
  camp[["label"]] <- Idents(camp)
  length(VariableFeatures(camp))
  
  anchors <- FindTransferAnchors(reference = camp, query = obj, 
                                 dims = seq(dims), reduction = "cca")
  predictions <- TransferData(anchorset = anchors, weight.reduction = "cca", 
                              refdata = camp$label, dims = seq(dims))
  
  obj <- AddMetaData(obj, metadata = predictions[,c(1,ncol(predictions))])
  
  return(obj)

}
