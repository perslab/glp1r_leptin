#' Create edgeR object
#'
#' @title
#' @import edgeR
#' @return
#' @author dylanmr
#' @export
edger_prep <- function(seur, trt_group = "geno", lib.size = 1e4,
                       filter_column = NULL, filter_value = NULL,
                       assay_name = "RNA", feature_pattern_exclude = "^Gm|Rik$", 
                       dataset = NULL, celltypes = celltypes, celltype_column = NULL,
                       min.count = 1, min.total.count = 5) {
  
  # Validate inputs
  stopifnot("Seurat" %in% class(seur), assay_name %in% names(seur@assays))
  DefaultAssay(seur) <- "RNA"
  
  # Subset cells based on metadata column and filter features
  cells_to_keep <- colnames(seur)[grepl(paste0("^",celltypes, "$"), seur[[celltype_column]][,1])]
  print(length(cells_to_keep))
  features_to_keep <- rownames(seur@assays[[assay_name]])[!grepl(feature_pattern_exclude, rownames(seur@assays[[assay_name]]))]
  seur <- subset(seur, cells = cells_to_keep, features = features_to_keep)
  
  if(!is.null(dataset)){
    cells_to_keep <- colnames(seur)[grepl(paste0("^",dataset, "$"), seur[["dataset"]][,1])]
    print(length(cells_to_keep))
    seur <- subset(seur, cells = cells_to_keep)
  }
  
  if(!is.null(filter_column)){
    cells_to_keep <- colnames(seur)[grepl(paste0("^",filter_value, "$"), seur[[filter_column]][,1])]
    print(length(cells_to_keep))
    seur <- subset(seur, cells = cells_to_keep)
  }
  
  y <- Seurat2PB(seur, sample = "hash.mcl.ID", cluster= trt_group)
  keep.samples <- y$samples$lib.size > lib.size
  y <- y[, keep.samples]
  meta <- seur[[]] %>% distinct(hash.mcl.ID, .keep_all = T) 
  meta <- meta[match(y$samples$sample, meta$hash.mcl.ID),]
  y$samples <- bind_cols(meta, y$samples)
  
  keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count = min.count, min.total.count = min.total.count)
  table(keep.genes)
  y <- y[keep.genes, , keep=FALSE]
  y <- normLibSizes(y)
  
  return(y)
  
}

#' wrapper to run scDist
#' @title
#' @param obj seurat object
#' @param fixed.effects character column with groups of interest (maximum two levels)
#' @param random.effects character vector of columns with batch information
#' @param assay seurat assay to use to calculate distance
#' @param ident_column character name of column of clusters
#' @param subset_cells regex format for clusters to keep
#' @param d numeric number of PC dimensions used to calculate distance
#' @param nfeats numeric number of genes used to calculate PCA
#' @return
#' @author dylanmr
#' @export
run_scdist <- function(obj, fixed.effects = "geno", assay="SCT", 
                            ident_column = "predicted.celltype", 
                            subset_cells = celltypes, 
                            random.effects = c("hash.mcl.ID", "orig.ident","sex_predicted","time", "treatment"),
                            d=20,
                            nfeats=5000) {

  DefaultAssay(obj) <- "RNA"
  obj <- process_seurat(obj, method = "qpoisson", cluster=F, nfeats = nfeats)
  cells <- colnames(obj)[grep(paste("^",subset_cells, "$",collapse="|", sep = ""), obj[[ident_column]][,1])]
  obj <- subset(obj, cells=cells)
  mat <- GetAssayData(obj, slot="scale.data", assay=assay) %>%  as.matrix()
  scd_res <- scDist(mat, meta.data=obj@meta.data, min.counts.per.cell = 5,
         fixed.effects = fixed.effects, 
         random.effects = random.effects, 
         clusters = ident_column, d = d)$results
  scd_res %<>% mutate(padj = p.adjust(p.val))
  return(scd_res)

}

#' Calculate DEG using edgeR
#'
#' @title
#' @param y pseudobulk obj generated from edger_prep
#' @param mm formula for differential expression
#' @param contrast vector to construct comparisons
#' @return
#' @author dylanmr
#' @export
run_edger <- function(y, mm = "~ 0 + cluster + sex_predicted + seq_pool", contrast = c(1,-1,0,0)) {

  design <- model.matrix(as.formula(mm), data = y$samples)
  colnames(design) <- gsub("/","_",colnames(design))
  
  y <- estimateDisp(y, design, robust=T)
  fit <- glmQLFit(y, design, robust=T)
  qlf <- glmQLFTest(fit, contrast=contrast)

}


#' Create Milo object
#' @title
#' @param obj seurat object
#' @param pca_reduction_name name of where PCA is stored in seurat object
#' @param umap_reduction_name name of where UMAP is stored in seurat object
#' @param k numeric number of neighbors
#' @param d numeric number of dimensions
#' @param prop 
#' @param refinement_scheme
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


