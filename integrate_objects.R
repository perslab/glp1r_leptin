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
          RunSexPrediction(., assay = "RNA") %>% 
          process_seurat(., method = "log", res = clustering_res, dims = clustering_dims, features = c("Tsix|Xist|Uty|Ddx3y"))
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
      ob_db_neuron_scenic,
      create_loom(ob_db_neuron, here::here("SCENIC/input/ob_db_neuron.loom"), subset_value=NULL)
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
      diet_neurons_scenic,
      create_loom(diet_neurons, here::here("SCENIC/input/diet_neurons.loom"), subset_value=NULL)
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
      glp1rleprko_neurons_scenic,
      create_loom(glp1rleprko_neurons, here::here("SCENIC/input/glp1rleprko_neurons.loom"), subset_value=NULL)
    ),
    tar_target(
      agrpleprko_neurons,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "agrpfl") %>% 
          filter_cells_by_column(., column_name = "l1_neur", value = "Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
    tar_target(
      agrpleprko_neurons_scenic,
      create_loom(agrpleprko_neurons, here::here("SCENIC/input/agrpleprko_neurons.loom"), subset_value=NULL)
    ),
    tar_target(
      campbell_scrna,
      qs::qread(here::here("external_data/20231027_ARC_filtered.qs"))
    ),
    tar_target(
      myers_scrna,
      prep_myers(here::here("external_data/myers_sctrap/"), min.features = 1000, min.cells=10)
    )
  )