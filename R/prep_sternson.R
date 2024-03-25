#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param path
#' @return
#' @author dylanmr
#' @export
prep_sternson <- function(path=NULL) {

  data <- data.table::fread(here::here("data/sternson/Agrp.Pomc.fedfast.counts.txt.gz")) %>% data.frame
  rownames(data) <- data$symbol
  meta <- str_split_fixed(colnames(data), "[.]", 3) %>% data.frame
  rownames(meta) <- colnames(data)
  dds <- DESeqDataSetFromMatrix(data, meta, ~X2)
  vsd <- vst(dds, blind = F)
  return(vsd)

}
