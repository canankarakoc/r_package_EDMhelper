#' ccm_summary
#'
#' Computes the summary statistics of the convergent cross mapping output.
#' @param ccm_output output of the convergent cross mapping.
#' @param predtype A character indicating which metric is used as predictive skill. Defaults to "rho". Can be "rho", "mae", or "rmse".
#' @keywords rEDM, ccm
#' @return list of mean, lower and upper confidence intervals of ccm iterations for the predictive skill used.
#' @import stats
#' @export

ccm_summary <- function (ccm_output, predtype="rho"){

  direction <- paste(ccm_output$target_column , "causes",  ccm_output$lib_column)
  ccm_output$predtype<-ccm_output[,predtype]

  summaryout   <-  with(ccm_output, aggregate(cbind(predtype=predtype),
                                       list(direction=direction,
                                            lib_size=lib_size),
                                       function(x) quantile(x, c(0.025,
                                                                 pnorm(-1,0,1), 0.5,
                                                                 pnorm(1,0,1), 0.975),
                                                            na.rm=T)))
  summaryout<-data.frame(summaryout[,1:2], unlist(summaryout[,3]))
  colnames(summaryout)[3:7]<-gsub("X", "Q", colnames(summaryout)[3:7], fixed=T)

  return(summaryout)
}


