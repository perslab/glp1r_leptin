#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param agrp
#' @param nameme1
#' @param nameme2
#' @return
#' @author dylanmr
#' @export
create_new_column_seurat <- function(obj, selected_columns = c("treatment", "time", "geno"), new_column_name = "group") {

  df <- obj[[]]
  
  # Check if all selected columns are in the dataframe
  if(!all(selected_columns %in% names(df))) {
    stop("Not all selected columns are in the dataframe")
  }
  
  # Create the new column by pasting the selected columns together
  # The 'apply' function is used to iterate over rows
  obj[[new_column_name]] <- apply(df[, selected_columns, drop = FALSE], 1, paste, collapse = "_")
  
  # Return the modified dataframe
  return(obj)

}
