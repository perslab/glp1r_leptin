# Seurat processing params
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

#targets options
tar_option_set(error="null")

analyze_lepip <- 
  list(
    tar_target(
      arc_ob_db,
      {
        classify_cells(ob_db_neuron, path = here::here("cell_markers/coarse_markers_vmh.txt"), column_name = "l2_arc_region", clustering_res = clustering_res,
                       clustering_dims = clustering_dims) %>% 
          filter_cells_by_column(., column_name = "l2_arc_region", value = "ARC_Neurons") %>% 
          process_seurat(., method = "integrate", batch = integration_batch,
                         res = clustering_res, dims = clustering_dims) %>% 
          map_ref(., dims=30, column_name = "internal_labs", ref = campbell_scrna)
      }
    ),
    tar_target(
      lepip_to_remove,
      union(names(which(prop.table(table(arc_ob_db$dataset, arc_ob_db$seurat_clusters), margin=2)[2,]>0.4)),
            names(which(tapply(arc_ob_db$prediction.score.max, arc_ob_db$seurat_clusters, mean)<0.7)))
    ),
    tar_target(
      arc_lepip,
      {
        subset(arc_ob_db, idents = lepip_to_remove, invert=T) %>% 
          filter_cells_by_column(., column_name = "internal_labs", 
                                 values = c("Lpar1_oligo","Hdc","Sim1/Rprm","Sim1/Ebf3"), invert = T) %>% 
          process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims, return_model=T) %>%
          create_new_column_seurat(., selected_columns = c("geno","treatment","time"), new_column_name = "group") 
      }
    ),
    tar_target(
      arc_lepip_labeled,
      arc_lepip %>% modify_column(., "seurat_clusters", c("0","13","26"),"0")
    ),
    tar_target(
      arc_lepip_h5, 
      convert_to_h5(arc_lepip_labeled, filename = here::here("h5_objects/arc_lepip.h5"))
    ),
    tar_target(
      arc_lepip_loom, 
      create_loom(arc_lepip_labeled, here::here("SCENIC/input/arc_lepip.loom"), subset_value=NULL)
    ),
    tar_target(
      celltypes,
      names(table(arc_lepip_labeled$seurat_clusters))[table(arc_lepip_labeled$seurat_clusters) > 1000]
    ),
    tar_target(
      pseudobulk_ob, 
      edger_prep(arc_lepip_labeled, dataset = "lepip", celltype_column="seurat_clusters", trt_group = "geno", celltypes = celltypes, filter_column = "treatment", filter_value = "Sal"),
      pattern = map(celltypes),
      iteration="list"
    ),
    tar_target(
      pseudobulk_lepip, 
      edger_prep(arc_lepip_labeled, dataset = "lepip", celltype_column="seurat_clusters", trt_group = "group", celltypes = celltypes),
      pattern = map(celltypes),
      iteration="list"
    ),
    tar_target(
      deg_ob,
      run_edger(pseudobulk_ob, mm = "~ 0 + cluster + seq_pool", contrast = c(1 , -1, 0)),
      pattern = map(pseudobulk_ob),
      iteration="list"
    ),
    tar_target(
      obwt_celldist, 
      {
        filter_cells_by_column(arc_lepip_labeled, column_name = "treatment", values = "Sal") %>% 
          subset(., subset = dataset == "lepip" & sex_predicted == 1 & seurat_clusters %in% celltypes) %>% 
          run_scdist(obj = ., fixed.effects = c("geno", "time", "seq_pool"),  
                     assay="SCT", ident_column = "seurat_clusters",
                     random.effects = c("orig.ident","hash.mcl.ID"), 
                     d = 20)
      }
    ),
    tar_target(
      lepip_celldist,
      {
        subset(arc_lepip_labeled, subset = time %in% c(1,3,6) & seurat_clusters %in% celltypes & dataset == "lepip") %>%
          create_new_column_seurat(., selected_columns = c("seurat_clusters","geno"), new_column_name = "scdist_group") %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "seq_pool","time"),
                     assay="SCT", ident_column = "scdist_group", 
                     random.effects = c("hash.mcl.ID","orig.ident"),
                     d = 20)
      }
    ),
    tar_target(
      focus_object,
      subset(arc_lepip_labeled, subset = seurat_clusters %in% c("0", "9", "26"))
    )
  )