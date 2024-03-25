#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
prep_dowsett <- function(method = "log") {

  hypomap_full <- readRDS(url("https://www.repository.cam.ac.uk/bitstream/handle/1810/340518/hypoMap.rds?sequence=2&isAllowed=y", "rb"))
  
  dowsett <- hypomap_full %>%
    subset(Batch_ID == 'Dowsett10xnuc_batch_1') %>%
    subset(Author_Class_Curated == "Neurons")
  
  #change default assay
  DefaultAssay(dowsett) <- 'RNA'
  
  #keep only RNA
  dowsett <- DietSeurat(dowsett, assays = "RNA") %>%
    #run SCT and PCA
    process_seurat(., method = method, res=0.8, dims = 30)

}
