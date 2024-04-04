# Seurat processing params
integration_batch <- "orig.ident"
clustering_res <- 0.8
clustering_dims <- 40

# targets options
tar_option_set(error = "null")

analyze_diet <-
  list(
    tar_target(
      arc_diet,
      {
        classify_cells(diet_neurons, path = here::here("cell_markers/coarse_markers_vmh.txt"), column_name = "l2_arc_region",
                       clustering_res = clustering_res, clustering_dims = clustering_dims) %>% 
        filter_cells_by_column(., column_name = "l2_arc_region", value = "ARC_Neurons") %>%
        process_seurat(., method = "integrate", batch = integration_batch, res = clustering_res, dims = clustering_dims) %>%
        project_umap(query = ., ref = arc_lepip_labeled, dims = 30, label_to_transfer="seurat_clusters", reference.assay = "integrated", query.assay = "integrated") 
      }
    ),
    tar_target(
      diet_to_remove,
      names(which(tapply(arc_diet$predicted.celltype.score, arc_diet$seurat_clusters, mean)<0.7))
    ),
    tar_target(
      arc_diet_labeled,
      subset(arc_diet, subset = seurat_clusters %in% diet_to_remove, invert=T)
    ),
    tar_target(
      celltypes_diet,
      names(table(arc_diet_labeled$predicted.celltype))[table(arc_diet_labeled$predicted.celltype) > 300]
    ),
    tar_target(
      hfdveh_celldist,
      {
        subset(arc_diet_labeled, subset = treatment %in% c("Chow", "HFD") & predicted.celltype %in% celltypes_diet) %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "hash_pool"), assay = "SCT", ident_column = "predicted.celltype", 
                   random.effects = c("orig.ident","hash.mcl.ID"), d = 20)
      }
    ),
    tar_target(
      fastveh_celldist,
      {
        subset(arc_diet_labeled, subset = treatment %in% c("Chow", "Fast") & predicted.celltype %in% celltypes_diet) %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "hash_pool"), assay = "SCT", ident_column = "predicted.celltype", 
                     random.effects = c("orig.ident","hash.mcl.ID"), d = 20)
      }
    ),
    tar_target(
      refeed_celldist,
      {
        subset(arc_diet_labeled, subset = treatment %in% c("Refeed", "Fast") & predicted.celltype %in% celltypes_diet) %>% 
          run_scdist(obj = ., fixed.effects = c("treatment", "hash_pool"), assay = "SCT", ident_column = "predicted.celltype", 
                     random.effects = c("orig.ident","hash.mcl.ID"), d = 20)
      }
    ),
    tar_target(
      pseudobulk_alldiet,
      edger_prep(arc_diet_labeled, celltype_column = "predicted.celltype", trt_group = "treatment", celltypes = celltypes_diet),
      pattern = map(celltypes_diet),
      iteration = "list"
    )
  )