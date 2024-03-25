# Seurat processing params
ann_column <- "ann_neuron"
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

knockout_agrp <- 
  list(
    tar_target(
      agrpko,
      integrated_neuron_agrpko %>%
        filter_cells_by_column(., column_name = "genetic", values = "Cre neg; Flox WT", invert=T) %>% 
        classify_cells_v2(., path = "/projects/dylan/target_testing/cell_markers/coarse_markers_vmh.txt") %>% 
        filter_cells_by_column(., column_name = "ann_neuron_v2", value = "ARC_Neurons") %>% 
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) %>% 
        map_ref(., dims=30, column_name = "internal_labs")
    ),
    tar_target(
      agrpko_scdist, 
      run_scdist(obj = agrpko, fixed.effects = "geno",  assay="SCT", column = "internal_labs",
                 random.effects = c("hash.mcl.ID", "hash_pool", "sex_predicted", "orig.ident"), 
                 d = 20, subset = celltypes)
    ),
    tar_target(
      celltypes_to_analyze,
      names(table(agrpko$newlabs))[table(glp1rko$newlabs)>200]
    ),
    tar_target(
      agrpko_pb,
      voom_wrapper(agrpko, filter_value = "Agrp", filter_column = "internal_labs", group_of_interest = "group", 
                   min_cells = 10, min_count_filter = 5, 
                   pseudobulk_var = c("hash.mcl.ID"))
    )
  )