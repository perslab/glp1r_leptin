#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
run_gex_qc <- function(obj, times, batch) {
  
  # convert to single cell experiment
  ds <- as.SingleCellExperiment(obj)
  ds <- scuttle::addPerCellQC(ds, percent_top=c(20,50,100,200))
  
  # calculate qc columns
  ds$total_features <- ds$detected
  ds$log10_total_features <- log10(ds$detected)
  ds$total_counts <- ds$sum
  ds$log10_total_counts <- log10(ds$sum+1)
  ds$featcount_ratio <- ds$log10_total_counts/ds$log10_total_features
  ds$featcount_dist <- getFeatCountDist(ds)
  ds$pct_counts_top_20_features <- colData(ds)[[intersect(c("percent_top_20","pct_counts_in_top_20_features","percent.top_20"), colnames(colData(ds)))[[1]]]]
  ds$pct_counts_top_50_features <- colData(ds)[[intersect(c("percent_top_50","pct_counts_in_top_50_features","percent.top_50"), colnames(colData(ds)))[[1]]]]
  
  # build matrix to track qc
  discard <- matrix(c(scuttle::isOutlier(ds$log10_total_counts, nmads = 5, type = "lower", batch = ds[[batch]]) , 
                      scuttle::isOutlier(ds$log10_total_counts,nmads = 2.5,type = "higher", batch = ds[[batch]]) ,
                      scuttle::isOutlier(ds$log10_total_features,nmads = 5,type = "lower", batch = ds[[batch]]) , 
                      scuttle::isOutlier(ds$log10_total_features, nmads = 2.5,type = "higher", batch = ds[[batch]]) ,
                      scuttle::isOutlier(ds$pct_counts_top_20_features, nmads = 5,type = "both", batch = ds[[batch]]) ,
                      scuttle::isOutlier(ds$pct_counts_top_50_features, nmads = 5,type = "both", batch = ds[[batch]]),
                      scuttle::isOutlier(ds$featcount_dist, nmads = 5,type = "both", batch = ds[[batch]]), 
                      ds$HTO_mcl_margin < (mean(ds$HTO_mcl_margin) - sd(ds$HTO_mcl_margin)*3)),
                    nrow = length(ds$log10_total_counts))

  # filter cells
  discard <- rowSums(discard)>=times
  cells_to_keep <- colnames(ds)[!discard]
  obj <- subset(obj, cells = cells_to_keep)

}



getFeatCountDist <- function(df, do.plot=FALSE, linear=TRUE){
  if(is(df,"SingleCellExperiment")) df <- as.data.frame(colData(df))
  if(linear){
    mod <- lm(log10_total_features~log10_total_counts, data=df)
    pred <- predict.lm(mod, newdata=data.frame(log10_total_counts=df$log10_total_counts))
  }else{
    mod <- loess(log10_total_features~log10_total_counts, data=df)
    pred <- predict.loess(mod, newdata=data.frame(log10_total_counts=df$log10_total_counts))
  }
  df$diff <- df$log10_total_features - pred
  if(do.plot){
    library(ggplot2)
    ggplot(df, aes(x=total_counts, y=total_features, colour=diff)) + 
      geom_point() + geom_smooth(method = "loess", col="black")
  }
  df$diff
}