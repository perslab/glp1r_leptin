# Seurat processing params
ann_column <- "ann_neuron"
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

integrate_objects <-
  list(
    tar_target(
      integrated_seurat,
      {
        merge_mats(filtered_seurat) %>% 
          process_seurat(., method = "log", res = clustering_res, dims = clustering_dims) %>% 
          RunSexPrediction(., assay = "RNA")
      }
    ),
    tar_target(
      ob_db_neuron,
      filter_cells_by_column(integrated_seurat, column_name = "dataset", value = c("reactivation","lepip")) %>% 
        filter_cells_by_column(., column_name = "geno", value = "react", invert = T) %>% 
        filter_cells_by_column(., column_name = "l1_neur", value = "Neurons") %>% 
        process_seurat(., method = "log", res = clustering_res, dims = clustering_dims)
    ),
    tar_target(
      ob_db_glia,
      filter_cells_by_column(integrated_seurat, column_name = "dataset", value = c("reactivation","lepip")) %>% 
        filter_cells_by_column(., column_name = "geno", value = "react", invert = T) %>% 
        filter_cells_by_column(., column_name = "l1_neur", value = "Neurons", invert = T) %>% 
        process_seurat(., method = "log", res = clustering_res, dims = clustering_dims)
    ),
    tar_target(
      diet_neurons,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "diet") %>% 
          filter_cells_by_column(., column_name = "l1_neur", value = "Neurons") %>% 
          process_seurat(., method = "log", res = clustering_res, dims = clustering_dims)
      }
    ),
    tar_target(
      diet_glia,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "diet") %>% 
          filter_cells_by_column(., column_name = "l1_neur", value = "Neurons", invert=T) %>% 
          process_seurat(., method = "log", res = clustering_res, dims = clustering_dims)
      }
    ),
    tar_target(
      glp1rleprko_neurons,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = c("glp1rleprko")) %>% 
          filter_cells_by_column(., column_name = "l1_neur", value = "Neurons") %>% 
          process_seurat(., method = "log", res = clustering_res, dims = clustering_dims)
      }
    ),
    tar_target(
      agrpleprko_neurons,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "agrpfl") %>% 
          filter_cells_by_column(., column_name = "l1_neur", value = "Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    )
  )



# 
# 
# 
# 
# 
# 
# 
# 