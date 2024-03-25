#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
run_scdist <- function(obj, fixed.effects = "geno", assay="SCT", ident_column = "predicted.celltype", 
                       subset_cells = celltypes, subset_column = "column", subset_groups = "groups",
                       random.effects = c("hash.mcl.ID", "orig.ident","sex_predicted","time", "treatment"), d=20, nfeats=5000) {

  DefaultAssay(obj) <- "RNA"
  obj <- process_seurat(obj, method = "qpoisson", cluster=F, nfeats = nfeats)
  cells <- colnames(obj)[grep(paste("^",subset_cells, "$",collapse="|", sep = ""), unlist(obj[[ident_column]]))]
  obj <- subset(obj, cells=cells)
  mat <- GetAssayData(obj, slot="scale.data", assay=assay) %>%  as.matrix()
  scd_res <- scDist(mat, meta.data=obj@meta.data, min.counts.per.cell = 5,
         fixed.effects = fixed.effects, 
         random.effects = random.effects, 
         clusters = ident_column, d = d)$results
  scd_res %<>% mutate(padj = p.adjust(p.val))
  return(scd_res)

}

