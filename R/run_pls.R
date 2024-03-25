#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
run_pls <- function(pb, ncomp = 2, features = 100) {
  
  vsd <- vst(pb, blind=F)
  nsamp <- table(colData(vsd)$final_labels, colData(vsd)$geno) %>%  data.frame() %>% pull(Freq)
  cut <- nsamp[which.min(nsamp)]-2
  new_df <- colData(vsd) %>% 
    data.frame() %>% 
    rownames_to_column("cell") %>% 
    group_by(geno, final_labels) %>% 
    slice_sample(n=cut)

  vsd <- vsd[,colnames(vsd)%in%new_df$cell]
  
  splsda <- mint.splsda(X = t(assay(vsd)), 
                                  Y = factor(vsd$geno), ncomp = ncomp,
                                  scale = T,
                                  study =  factor(vsd$final_labels), 
                                  keepX = c(features)) 
  return(splsda)
}
