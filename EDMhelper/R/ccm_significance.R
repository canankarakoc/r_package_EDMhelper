#' significance_test
#'
#' Significance test for convergent cross mapping
#' @param ccm_summary output of a call from "ccm_summary"
#' @param predtype A character indicating which metric to be used for determining predictive skill. Defaults to "rho". Can be "rho", "mae", or "rmse".
#' @keywords rEDM, significance
#' @return
#' @import
#' @export

significance_test <- function (ccm_summary, predtype="rho") {

  significance_data <- as.data.frame(ccm_summary$direction)


  significance_data$significance <- ifelse(ccm_summary[ccm_summary$lib_size==max(ccm_summary$lib_size),]$rho[,1]>0 &
                      (ccm_summary[ccm_summary$lib_size==max(ccm_summary$lib_size),]$predtype[,1]-
                         ccm_summary[ccm_summary$lib_size==min(ccm_summary$lib_size),]$predtype[,1])>0, "1","0")

return(significance_data)
}

