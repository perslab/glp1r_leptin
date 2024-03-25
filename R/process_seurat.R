#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param .
#' @param method
#' @param res
#' @param features
#' @param dims
#' @return
#' @author dylanmr
#' @export

process_seurat <- function(obj, method, ref_datasets = NULL, k.anchor=5,k.weight=100,
                           res=NULL, features=NULL, dims=NULL, 
                           batch=NULL, return_model = FALSE, cluster=TRUE, 
                           type = "seur", nfeats = 2000, neighbor=FALSE) {
  
  if (is(obj, "SingleCellExperiment")) {
    obj <- CreateSeuratObject(counts = counts(obj))
  }
  
  if(!is.null(features)) {
    features_to_remove <- rownames(obj)[grepl(features, rownames(obj))]
    obj <- subset(obj, features = setdiff(rownames(obj), features_to_remove))
  }
  
  obj <- switch(
    method,
    integrate = process_integrate(obj, batch, nfeats, ref_datasets, k.anchor, k.weight),
    log = process_log(obj, nfeats, batch),
    glm = process_glm(obj, nfeats, batch),
    qpoisson = process_qpoisson(obj, nfeats, batch),
    stop(paste("Unknown method:", method))
  )
  
  if (cluster) {
    obj <- cluster_obj(obj, dims, return_model, neighbor, res)
  }
  
  return(obj)
}

process_integrate <- function(obj, batch, nfeats, ref_datasets, k.anchor, k.weight) {
  if(is.list(obj)) {
    obj_list <- obj
  } else {
    DefaultAssay(obj) <- "RNA"
    obj_list <- SplitObject(obj, split.by = batch)
  }
  
  obj_list <- lapply(X = obj_list, FUN = function(x) {
    if(sum("SCT" %in% names(x@assays))>0) {
      x[["SCT"]] <- NULL
    }
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeats)
  })
  
  features <- SelectIntegrationFeatures(object.list = obj_list)
  
  obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  anchors <- FindIntegrationAnchors(object.list = obj_list, reference = ref_datasets, k.anchor = k.anchor, anchor.features = features, reduction = "rpca")
  integrated <- IntegrateData(anchorset = anchors, k.weight=k.weight)
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated) %>% RunPCA(.)
  return(integrated)
  
}

process_log <- function(obj, nfeats, batch) {
  if (inherits(obj, "list")) {
    obj <- obj[[1]]
  }
  
  DefaultAssay(obj) <- "RNA"
  
  NormalizeData(obj) %>%
    FindVariableFeatures(., selection.method = "vst", nfeatures = nfeats) %>%
    ScaleData(., vars.to.regress=batch) %>%
    RunPCA(.)
}

process_glm <- function(obj, nfeats, batch) {
  if (inherits(obj, "list")) {
    obj <- obj[[1]]
  }
  
  SCTransform(obj, method = "glmGamPoi", batch_var=batch, variable.features.n = nfeats) %>%
    RunPCA(.)
}

process_qpoisson <- function(obj, nfeats, batch) {
  if (inherits(obj, "list")) {
    obj <- obj[[1]]
  }
  
  SCTransform(obj, method = "qpoisson", variable.features.n = nfeats, vars.to.regress = batch) %>%
    RunPCA(.)
}

cluster_obj <- function(obj, dims, return_model, neighbor, res) {
  obj %>%
    RunUMAP(., dims = seq(dims), return.model=return_model) %>%
    FindNeighbors(., dims = seq(dims), return.neighbor=neighbor) %>%
    FindClusters(., resolution = res)
}