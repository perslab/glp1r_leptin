#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param seur
#' @param mincells
#' @return
#' @author dylanmr
#' @export
#' 
#' 

voom_wrapper <- function(seurat_obj, 
                         assay_name = "RNA", 
                         feature_pattern_exclude = "^Gm|Rik$", 
                         min_cells = 5, min_reps = 2,
                         filter_column = NULL, 
                         filter_value = NULL,
                         pseudobulk_var = c("hash.mcl.ID"), 
                         min_count_filter = 5,
                         group_of_interest = NULL) {
  
  # Validate inputs
  stopifnot("Seurat" %in% class(seurat_obj), assay_name %in% names(seurat_obj@assays))
  
  # Subset cells based on metadata column and filter features
  cells_to_keep <- colnames(seurat_obj)[grepl(filter_value, seurat_obj[[filter_column]][,1])]
  features_to_keep <- rownames(seurat_obj@assays[[assay_name]])[!grepl(feature_pattern_exclude, rownames(seurat_obj@assays[[assay_name]]))]
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep, features = features_to_keep)
  
  # Process grouping columns for pseudobulk
  if(length(pseudobulk_var>1)) {
    seurat_obj@meta.data %<>%  unite(., "pb_var", {{pseudobulk_var}}, remove = F)
  } else {
    seurat_obj@meta.data$pb_var <- seurat_obj@meta.data[{{pseudobulk_var}}]
  }
  
  # Unify metadata group_of_interest
  if(length(group_of_interest>1)) {
    seurat_obj@meta.data %<>%  unite(., "group_of_interest", {{group_of_interest}}, remove = F)
  } else {
    seurat_obj@meta.data$group_of_interest <- seurat_obj@meta.data[{{group_of_interest}}]
  }
  
  # Filter groups by cell count
  pb_sizes <- table(seurat_obj$pb_var)
  samples_keep <- names(pb_sizes[pb_sizes > min_cells])
  seur_filtered_ncells <- subset(seurat_obj, subset = pb_var %in% samples_keep)
  
  # Create a table counting the number of unique animals in each group
  samps_per_group <- table(seur_filtered_ncells$pb_var, seur_filtered_ncells$group_of_interest)
  groups_keep <- colnames(samps_per_group)[apply(samps_per_group, 2, function(x) sum(x > 0)) >= min_reps]
  seur_filtered_complete <- subset(seur_filtered_ncells, subset = group_of_interest %in% groups_keep)
  
  # Calculate pseudobulk matrix
  pseudobulk_matrix <- create_pseudobulk_matrix(seur_filtered_complete, assay_name)
  
  # Create metadata for voom
  metadata_for_voom <- distinct(seur_filtered_complete@meta.data, pb_var, .keep_all = T)
  metadata_for_voom <- metadata_for_voom[match(colnames(pseudobulk_matrix),metadata_for_voom$pb_var),]
  metadata_for_voom <- inner_join(metadata_for_voom, enframe(pb_sizes, name = "pb_var", value = "ncells") %>% mutate(ncells = as.numeric(ncells)), by = "pb_var")
  
  # Create a DGEList object
  y <- edgeR::DGEList(counts = pseudobulk_matrix, group = metadata_for_voom$group_of_interest, samples = metadata_for_voom)
  
  # Filter the expressions
  keep.exprs <- edgeR::filterByExpr(y, group = y$samples$group_of_interest, min.count = min_count_filter)
  print(table(keep.exprs))
  
  # Subset the DGEList object based on filtered expressions
  y <- y[keep.exprs,, keep.lib.sizes = FALSE]
  
  # Calculate normalization factors
  y <- edgeR::calcNormFactors(y, method="TMM")

  return(y)
}

# Helper functions that handle specific tasks within voom_wrapper

create_pseudobulk_matrix <- function(seurat_obj, assay = "RNA") {
  stopifnot("Seurat" %in% class(seurat_obj), assay %in% names(seurat_obj@assays))
  samples <- unique(seurat_obj$pb_var)
  # Retrieve counts from the specified assay
  assay_counts <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
  
  # Check if pb_var column exists in metadata
  if(!"pb_var" %in% colnames(seurat_obj@meta.data)) {
    stop("pb_var column does not exist in the Seurat object's metadata.")
  }
  
  # Calculate row sums for each group
  pbmat <- sapply(samples, function(group) {
    group_cells <- colnames(seurat_obj)[seurat_obj$pb_var == group]
    if(length(group_cells) == 0) {
      stop(paste("Group", group, "has no cells in the Seurat object."))
    }
    group_counts <- assay_counts[, group_cells, drop = FALSE]
    rowSums(group_counts)
  })
  
  return(pbmat)
}

