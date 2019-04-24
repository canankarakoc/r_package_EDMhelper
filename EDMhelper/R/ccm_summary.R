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

  summary   <-  with(ccm_output, aggregate(cbind(predtype=predtype),
                                       list(direction=direction,
                                            lib_size=lib_size),
                                       function(x) quantile(x, c(0.025,
                                                                 pnorm(-1,0,1), 0.5,
                                                                 pnorm(1,0,1), 0.975),
                                                            na.rm=T)))


}


