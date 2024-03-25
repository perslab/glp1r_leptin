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
      integrated_all_lepip,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "lepip") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
    tar_target(
      all_lepip_milo,
      create_milo(integrated_all_lepip)
    ),
    tar_target(
      all_neurons,
      {
        filter_cells_by_column(integrated_seurat, column_name = "ann_neuron", value = "Neurons") %>%
          filter_cells_by_column(., column_name = "dataset", value = "glp1rsun", invert=T) %>% 
          process_seurat(., method = "log", res = clustering_res, dims = clustering_dims) 
      }
    ),
    tar_target(
      integrated_neuron_lepip,
      {
        filter_cells_by_column(all_neurons, column_name = "dataset", value = "lepip") %>% 
        filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
    tar_target(
      integrated_arc_lepip,
      {
        classify_cells_v2(integrated_neuron_lepip, path = "/projects/dylan/target_testing/cell_markers/coarse_markers_vmh.txt") %>% 
          filter_cells_by_column(., column_name = "ann_neuron_v2", value = "ARC_Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) %>%
          create_new_column_seurat(., selected_columns = c("geno","treatment","time"), new_column_name = "group") %>% 
          map_ref(., dims=30, column_name = "internal_labs")
      }
    ),
    tar_target(
      integrated_glia_lepip,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "lepip") %>%
          filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons", invert = T) %>%
          create_new_column_seurat(., selected_columns = c("geno","treatment","time"), new_column_name = "group") %>% 
          process_seurat(., method = "integrate", batch = "hash_pool", res = clustering_res, dims = 20)
      }
    ),
    tar_target(
      integrated_neuron_agrpko,
      {
        filter_cells_by_column(all_neurons, column_name = "dataset", value = "agrpfl") %>% 
          filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
    tar_target(
      integrated_neuron_glp1rko,
      {
        filter_cells_by_column(all_neurons, column_name = "dataset", value = "glp1rfl") %>% 
          filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
    tar_target(
      integrated_neuron_diet,
      {
        filter_cells_by_column(all_neurons, column_name = "dataset", value = "diet") %>% 
          filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
    tar_target(
      integrated_arc_diet,
      {
        classify_cells_v2(integrated_neuron_diet, path = "/projects/dylan/target_testing/cell_markers/coarse_markers_vmh.txt") %>% 
          filter_cells_by_column(., column_name = "ann_neuron_v2", value = "ARC_Neurons") %>% 
          map_ref(., dims=30, column_name = "internal_labs") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims)
      }
  ),
    tar_target(
      integrated_glia_diet,
      {
        filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "diet") %>% 
          filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons", invert=T) %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
    # tar_target(
    #   integrated_neuron_dev,
    #   {
    #     filter_cells_by_column(all_neurons, column_name = "dataset", value = "dev") %>% 
    #       filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
    #       process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims, ref_datasets = sample(1:35, 20)) 
    #   }
    # ),
    tar_target(
      integrated_neuron_lepan,
      {
        filter_cells_by_column(all_neurons, column_name = "dataset", value = "lepan") %>% 
          filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
      }
    ),
  tar_target(
    integrated_glp1rsun,
    {
      filter_cells_by_column(integrated_seurat, column_name = "dataset", value = "glp1rsun") %>% 
        filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) 
    }
  ),
    tar_target(
      arc_ob_db,
      {
        filter_cells_by_column(all_neurons, column_name = "dataset", value = c("reactivation","lepip")) %>% 
          filter_cells_by_column(., column_name = "ann_neuron", value = "Neurons") %>% 
          filter_cells_by_column(., column_name = "geno", value = "React", invert = T) %>% 
          classify_cells_v2(integrated_neuron_lepip, path = "/projects/dylan/target_testing/cell_markers/coarse_markers_vmh.txt") %>% 
          filter_cells_by_column(., column_name = "ann_neuron_v2", value = "ARC_Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch, features = c("Tsix|Xist|Uty|Ddx3y"),
                         res = clustering_res, dims = clustering_dims) %>%
          map_ref(., dims=30, column_name = "internal_labs")
      }
    )
  )
  