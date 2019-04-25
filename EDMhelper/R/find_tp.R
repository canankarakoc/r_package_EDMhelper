#' find_tp
#'
#' Finds the prediction lag tp corresponding to a convergent cross mapping significance test.
#' @param sgn output of a call from "ccm_significance".
#' @param cutoff cuttoff for p-value. Default is 0.05.
#' @keywords rEDM, prediction lag
#' @return Data frame with the tps corresponding to the best predictive skill.
#' @export

find_tp<-function(sgn, cutoff=0.05) {

  best_tp <- NULL
  dirlst<-sort(unique(sgn$direction))

  for (i in 1:length(dirlst)){
    sbs<-which(sgn$pval[sgn$direction==dirlst[i]]<=cutoff)
    if(length(sbs)==0) {
      sbs<-which(sgn$tp[sgn$direction==dirlst[i]]==0)
    } else {
      sbs<-sbs[which.max(sgn$meanmaxL[sgn$direction==dirlst[i]][sbs])]

    }
    best_tp<-rbind(best_tp, sgn[sgn$direction==dirlst[i],][sbs,])
  }

  return(best_tp)
}





