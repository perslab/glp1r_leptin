tar_option_set(error="null")

analyze_lepip <- 
  list(
    tar_target(
      arc_ob_db_final, 
      subset(arc_ob_db, idents = c(35, 14, 45, 42, 33, 18, 19), invert=T) %>% 
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) %>% 
        modify_column(., "geno", "Crectrl", "WT") %>% 
        modify_column(., "geno", "KO", "db/db") %>% 
        modify_column(., "geno", "ko", "ob/ob") %>% 
        modify_column(., "dataset", "reactivation", "db/db") %>% 
        modify_column(., "dataset", "lepip", "ob/ob") 
    ),
    tar_target(
      arc_ob_db_notrt, 
      subset(arc_ob_db_final, subset = treatment != "Lep") %>% 
        process_seurat(., method = "integrate", batch = integration_batch, 
                       res = clustering_res, dims = clustering_dims)
    ),
    tar_target(
      arc_ob_db_notrt_int, 
      process_seurat(arc_ob_db_notrt, method = "integrate", batch = "geno", 
                     res = clustering_res, dims = clustering_dims) %>% 
        create_new_column_seurat(., c("orig.ident", "hash.mcl.ID"), "re")
    ),
    tar_target(
      arcobdb_h5, 
      convert_to_h5(arc_ob_db_notrt_int, filename = "/projects/dylan/target_testing/python_objects/arc_ob_db_notrt_int.h5")
    ),
    tar_target(
      celltypes,
      names(table(arc_ob_db_notrt_int$seurat_clusters))[table(arc_ob_db_notrt_int$seurat_clusters) > 600]
    ),
    tar_target(
      pseudobulk_ob, 
      edger_prep(arc_ob_db_notrt_int, dataset = "ob/ob", celltype_column="seurat_clusters", celltypes = celltypes),
      pattern = map(celltypes),
      iteration="list"
    ),
    tar_target(
      deg_ob,
      run_edger(pseudobulk_ob, mm = "~ 0 + cluster + sex_predicted + seq_pool", contrast = c(1 , -1, 0, 0)),
      pattern = map(pseudobulk_ob),
      iteration="list"
    ),
    tar_target(
      pseudobulk_db, 
      edger_prep(arc_ob_db_notrt_int, dataset = "db/db", celltype_column="seurat_clusters", celltypes = celltypes),
      pattern = map(celltypes),
      iteration="list"
    ),
    tar_target(
      deg_db,
      run_edger(pseudobulk_db, mm = "~ 0 + cluster + sex_predicted + hash_pool", contrast = c(1 , -1, 0, 0)),
      pattern = map(pseudobulk_db),
      iteration="list"
    ),
    tar_target(
      pseudobulk_combined, 
      arc_ob_db_notrt_int %>% 
        add_new_annotation(obj = .,new_name = "de_label", name = "geno", string = "b") %>% 
        edger_prep(., dataset = NULL, trt_group = "de_label", 
                   celltype_column="seurat_clusters", celltypes = celltypes),
      pattern = map(celltypes),
      iteration="list"
    ),
    tar_target(
      deg_combined,
      run_edger(pseudobulk_combined, mm = "~ 0 + cluster + sex_predicted + seq_pool", contrast = c(1 , -1, 0, 0, 0)),
      pattern = map(pseudobulk_combined),
      iteration="list"
    ),
    tar_target(
      obwt_celldist, 
      {
        run_scdist(obj = subset(arc_ob_db_notrt_int, subset = dataset == "ob/ob" & sex_predicted == 1), 
                     fixed.effects = c("geno", "time", "seq_pool"),
                     assay="SCT", ident_column = "seurat_clusters",
                     random.effects = c("re"), 
                     d = 20, subset_cells = celltypes)
      }
    ),
    tar_target(
      dbwt_celldist, 
      {
          run_scdist(obj =  subset(arc_ob_db_notrt_int, subset = dataset == "db/db"), 
                     fixed.effects = c("geno", "sex", "hash_pool"),
                     assay="SCT", ident_column = "seurat_clusters",
                     random.effects = c("re"), 
                     d = 20, subset_cells = celltypes)
      }
    ),
    tar_target(
      arc_ob_lepip, 
      subset(arc_ob_db_final, subset = dataset == "ob/ob") %>% 
        process_seurat(., method = "integrate", batch = integration_batch, 
                       res = clustering_res, dims = clustering_dims)
    ),
    tar_target(
      celltypes_lepip,
      names(table(arc_ob_lepip$seurat_clusters))[table(arc_ob_lepip$seurat_clusters) > 600]
    ),
    tar_target(
      ob_lepip_celldist,
        run_scdist(obj = subset(arc_ob_lepip, subset = dataset == "ob/ob" & time!=24 & geno == "ob"),
                   fixed.effects = c("treatment", "time", "seq_pool"),
                   assay="SCT", ident_column = "seurat_clusters",
                   random.effects = c("hash.mcl.ID"),
                   d = 20, subset_cells = celltypes_lepip)
    ),
    tar_target(
      wt_lepip_celldist,
      run_scdist(obj = subset(arc_ob_lepip, subset = dataset == "ob/ob" & time!=24 & geno == "wt"),
                 fixed.effects = c("treatment", "time", "seq_pool"),
                 assay="SCT", ident_column = "seurat_clusters",
                 random.effects = c("hash.mcl.ID"),
                 d = 20, subset_cells = celltypes_lepip)
    ),
    tar_target(
      focus_object,
      subset(arc_ob_db_notrt, subset = internal_labs %in% c("Agrp", "Trh/Cxcl12", "Pomc/Anxa2"))
    )
    # tar_target(
    #   lepip_celldist_ob, 
    #   subset(lepip, subset = geno == "ob" & time != 24) %>% 
    #     run_scdist(obj = ., fixed.effects = "treatment",  assay="SCT", column = "newlabs",
    #                random.effects = c("hash.mcl.ID", "seq_pool", "time", "orig.ident"), 
    #                d = 20, subset = celltypes)
    # ),
    # tarchetypes::tar_map(
    #   list(method = c("pca", "IEG")),
    #   tar_target(
    #     lepip_activity,
    #     predict_activity(method = method, lepip, ref_data = hipp_activity, ngene=100)
    #   )
    # )
  )