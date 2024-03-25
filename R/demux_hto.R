#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export

demux_hto <- function(obj, outs = "sc_preprocessing/hto_qc/", assay = "HTO") {
  outs <- here::here(outs)
  lane <- as.character(unique(obj$orig.ident))
  obj <- NormalizeData(obj, assay = assay, normalization.method = 'CLR', margin = 2, verbose = F)
  # Setting the threshold based on the quantiles of the negative and the positive clusters.
  q_l = 1
  q_h = 0.001
  obj <- HTODemux.mcl(obj, q_l =  q_l, q_h = q_h, assay=assay)
  table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Doublet']
  cat('#\t-\t\tDoublets:\t\t',table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Doublet'],'\n',
      '#\t-\t\tNegatives:\t\t',table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Negative'],'\n',
      '#\t-\t\tSinglets:\t\t',table(obj$HTO_mcl_classification.global)[names(table(obj$HTO_mcl_classification.global)) == 'Singlet'],'\n',
      sep = '', file = paste0(outs, lane,".log"))
  HTODemux_mcl.visualization(obj, q_l =  q_l, q_h = q_h, outs = outs, lane=lane)
  return(obj)
  
}

is_multimodal <- function(x, p.cutoff = 1e-2) {
  # Test if the expression distribution violates unimodal distribution.
  p = diptest::dip.test(x)$p.value
  return(p < p.cutoff)
}

select_hash_cutoff_mcl <- function(x, q_l = 1, q_h = 0.01, seed = 42) {
  # Model HTO data as a mixture of two Gaussian distributions (for normalized [across cells] data)
  # And select HTO cutoff based on mclust (Model based clustering).
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  # q_l and q_h are the quantiles for negative and postive cluster, respectively. 
  cutoff <- max(quantile(x[cl == low_cl], q_l), quantile(x[cl == high_cl], q_h))
  # The higher the cut off, the less false positive (the more false negative).
  return(cutoff)
}

hash_mcl_p <- function(x, seed = 3030, q_l = 1, q_h = 0.001) {
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  p.high_cl <- km$z[,high_cl]
  # Correct assignment error from Mclust
  p.high_cl[which(x < max(quantile(cl_center[low_cl], q_l), quantile(cl_center[low_cl], q_h)))] = 0
  names(p.high_cl) = names(x)
  return(p.high_cl)
}

HTO_classifcation = function(discrete, hto_mcl.p, assay){
  # Based on HTODemux (Seurat)
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = discrete)
  hash.max <- apply(X = hto_mcl.p, MARGIN = 2, FUN = max) # This returns the probability of the most likely HashID (based on the Hashtag distribution among cells)
  hash.maxID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = which.max)])
  hash.second <- apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) sort(x,decreasing = T)[2])
  hash.secondID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) order(x,decreasing = T)[2])])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, 'mcl', c("maxID", "secondID", "margin", "classification", "classification.global"), sep = "_")
  return(classification.metadata)
}

HTODemux.mcl <- function(object, assay = "HTO", q_l = 1, q_h = 0.005, seed = 42){
  # A function to find the threshold for each hastag, the P, and singlet, doublets and negative.
  # The input is the HTO data matrix (normalized across cells).
  
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x, q_l = q_l, q_h = q_h), future.seed = T))
  hto_mcl.cutoff$Multi_modal = apply(data,1,function(x) is_multimodal(x))
  print(hto_mcl.cutoff)
  hto_mcl.p = t(apply(data,1,function(x) hash_mcl_p(x, seed = seed, q_l = q_l, q_h = q_h)))
  discrete <- data
  discrete[discrete > 0] <- 0
  for (iter in rownames(x = data)) {
    values <- data[iter, ]
    cutoff <- hto_mcl.cutoff[iter,'cut_off']
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
  }
  classification.metadata <- HTO_classifcation(discrete, hto_mcl.p, assay)
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste(assay, "mcl", "classification", sep = '_')
  doublets <- rownames(x = object[[]])[which(object[[paste(assay, "mcl","classification.global", sep = "_")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  Idents(object) = factor(Idents(object), levels = c('Doublet', 'Negative', rownames(object@assays$HTO)))
  object$hash.mcl.ID <- Idents(object = object)
  return(object)
}

HTODemux_mcl.visualization <- function(object, assay = "HTO", q_l = 1, q_h = 0.005, seed = 42, outs, lane){
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto.data.wide = data.frame(t(data))
  hto.data.long = data.table::melt(data.table::setDT(hto.data.wide,keep.rownames = T), id.vars = 'rn', value.name = 'Expression',variable.name = 'hto')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x, q_l = q_l, q_h = q_h), future.seed = T), hto = colnames(hto.data.wide)[-1])
  
  p <- ggplot(hto.data.long, aes(x = Expression)) +
    geom_histogram(bins = 100) +
    geom_vline(data = hto_mcl.cutoff, aes(xintercept = cut_off), col = 'red') +
    facet_wrap(~hto,scales = 'free',ncol = 2) +
    xlab('Expression') +
    ylab('Counts') +
    scale_y_sqrt() +
    ggtitle('Individual HTO distributions', subtitle = paste0('Q_negative = ', q_l, '; Q_positive = ', q_h)) +
    theme_minimal(base_size = 10)
  
  ggsave(filename = paste0(lane,'_HTO-hist.png'),
         plot = p,
         path = outs,
         width = 10,
         height = 12)
}
