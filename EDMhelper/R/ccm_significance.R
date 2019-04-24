#' significance_test
#'
#' Significance test for convergent cross mapping
#' @param ccm_output output of a call from "ccm_easy"
#' @param predtype A character indicating which metric to be used for determining predictive skill. Defaults to "rho". Can be "rho", "mae", or "rmse".
#' @keywords rEDM, significance
#' @return Results from significance test.
#' @export

significance_test <- function (ccm_output, predtype="rho") {

  ccm_output$direction<-paste(ccm_output$target_column , "causes",  ccm_output$lib_column)

  significance_data <- ccm_output$direction
  significance_data_lst<-sort(unique(significance_data))
  pval<-data.frame(direction=significance_data_lst, pval=NA)

  for(i in 1:length(significance_data_lst)) {
    sbs<-which(significance_data==significance_data_lst[i])
    mx<-ccm_output[sbs,][ccm_output$lib_size[sbs]==max(ccm_output$lib_size[sbs]),predtype]
    mn<-ccm_output[sbs,][ccm_output$lib_size[sbs]==min(ccm_output$lib_size[sbs]),predtype]

    if(predtype=="rho") {
      pval[i,2]<-mean(mx<0 | (mx<=mn))
    } else {
      pval[i,2]<-mean(mx>=mn)
    }
  }

  return(pval)
}

