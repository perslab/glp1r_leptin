#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
#'

filter_cells_by_column <- function(input, column_name, values, invert = FALSE) {
  
  if (!is.vector(values) || is.null(values)) {
    stop("The 'values' parameter should be a non-null vector.")
  }
  
  if (inherits(input, "Seurat")) {
    seurat_list <- list(input)
  } else if (is.list(input) && all(sapply(input, inherits, "Seurat"))) {
    seurat_list <- input
  } else {
    stop("Input must be a Seurat object or a list of Seurat objects.")
  }
  
  filtered_objects <- lapply(seurat_list, function(s) {
    # Get cells that match any of the desired values in the specified column
    if (invert) {
      cells_to_keep <- !s@meta.data[[column_name]] %in% values
    } else {
      cells_to_keep <- s@meta.data[[column_name]] %in% values
    }
    
    # Subset the Seurat object to keep only the desired cells
    s_subset <- subset(s, cells = colnames(s)[cells_to_keep])
    
    return(s_subset)
  })
  
  # Check the length of the list and return appropriately
  if (length(filtered_objects) == 1) {
    return(filtered_objects[[1]])
  } else {
    return(filtered_objects)
  }
  
  return(filtered_objects)
  
}

