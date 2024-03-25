#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export

create_milo <- function(obj, pca_reduction_name = 'pca', umap_reduction_name = 'umap',
                            k = 60, d = 40, prop = 0.05, refinement_scheme = "graph") {
 
   # Convert Seurat object to SingleCellExperiment object
  sce <- as.SingleCellExperiment(obj)
  
  # Add PCA and UMAP embeddings
  reducedDim(sce, pca_reduction_name, withDimnames = TRUE) <- obj[[pca_reduction_name]]@cell.embeddings
  reducedDim(sce, umap_reduction_name, withDimnames = TRUE) <- obj[[umap_reduction_name]]@cell.embeddings
  
  # Create Milo object
  milo.obj <- Milo(sce)
  milo.obj <- buildGraph(milo.obj, k = k, d = d)
  milo.obj <- makeNhoods(milo.obj, prop = prop, k = k, d = d, refined = TRUE, refinement_scheme = refinement_scheme)
  
  # Add metadata and count cells
  colData(milo.obj)$ObsID <- paste(colData(milo.obj)$orig.ident, colData(milo.obj)$hash.mcl.ID, sep = "_")
  milo.obj <- countCells(milo.obj, meta.data = data.frame(colData(milo.obj)), samples = "ObsID")
  
  # Return milo.obj
  return(milo.obj)
}
