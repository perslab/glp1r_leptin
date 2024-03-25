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
map_ref <- function(obj = obj, dims = 30, 
                    ref = "/projects/amj/my_projects/reference_ARCDVC/20231010_full_integration/1_dat/1_seurat_obj/20231027_ARC_filtered.qs", 
                    column_name = "final_labels") {

  ref_obj <- qs::qread(ref)
  ref_sub <- subset(ref_obj, cells = sample(Cells(ref_obj), 20000))
  ref_sub$stablecamp7 <- ifelse(grepl("Agrp",ref_sub$stablecamp7), "Agrp", ref_sub$stablecamp7)
  anchors <- FindTransferAnchors(reference = ref_sub, query = obj, dims = seq(dims), reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = ref_sub$stablecamp7, dims = seq(dims))
  obj[[column_name]] <- predictions$predicted.id
  obj[[column_name]] <- ifelse(grepl("Agrp", obj[[column_name]][,1]), "Agrp", obj[[column_name]][,1])
  obj[["prediction.score.max"]] <- predictions$prediction.score.max
  return(obj)

}


