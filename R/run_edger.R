#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
run_edger <- function(y, mm = "~ 0 + cluster + sex_predicted + seq_pool", contrast = c(1,-1,0,0)) {

  design <- model.matrix(as.formula(mm), data = y$samples)
  colnames(design) <- gsub("/","_",colnames(design))
  
  y <- estimateDisp(y, design, robust=T)
  fit <- glmQLFit(y, design, robust=T)
  qlf <- glmQLFTest(fit, contrast=contrast)

}
