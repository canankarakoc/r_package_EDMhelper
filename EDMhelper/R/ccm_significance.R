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
  tpdat<-ccm_output$tp
  significance_data_lst<-sort(unique(significance_data))
  tp_lst<-sort(unique(ccm_output$tp))

  pval<-data.frame(direction=rep(significance_data_lst, each=length(tp_lst)), tp=rep(tp_lst, length(significance_data_lst)), pval=NA, meanmaxL=NA)

  n<-1
  for(i in 1:length(significance_data_lst)) {
    for(j in 1:length(tp_lst)) {
      sbs<-which(significance_data==significance_data_lst[i] & tpdat==tp_lst[j])
      mx<-ccm_output[sbs,][ccm_output$lib_size[sbs]==max(ccm_output$lib_size[sbs]),predtype]
      mn<-ccm_output[sbs,][ccm_output$lib_size[sbs]==min(ccm_output$lib_size[sbs]),predtype]

      if(predtype=="rho") {
        pval[n,3]<-mean(mx<0 | (mx<=mn), na.rm=T)
      } else {
        pval[n,3]<-mean(mx>=mn, na.rm=T)
      }
      pval[n,4]<-mean(mx, na.rm=T)
      n<-n+1
    }
  }

  return(pval)
}

