#' find_tp
#'
#' Finds the prediction lag tp corresponding to a convergent cross mapping significance test.
#' @param pval output of a call from "ccm_significance".
#' @param cutoff cuttoff for p-value. Default is 0.05.
#' @keywords rEDM, prediction lag
#' @return Data frame with the tps corresponding to the best predictive skill.


find_tp<-function(sgn, cutoff=0.05) {

  best_tp <- list()

  for (i in levels(sgn$direction)){

    for (j in length(sgn$pval)){

      if(sgn$pval[j] <= cutoff) {

        best_tp[[i]]=sgn$tp[which.min(sgn$pval)]

      } else {

        best_tp[[i]] = 0

      }
    }
  }

  return(best_tp)
}





