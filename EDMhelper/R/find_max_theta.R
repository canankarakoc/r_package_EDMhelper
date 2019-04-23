#' find_max_theta
#'
#' Finds the embedding dimension corresponding to the best predictive skill for a simplex call.
#' @param smap_out output of a call from "smap"
#' @param buffer A number from 0 to 1, indicating the allowed buffer distance between the "best" E, and the E that is used. If mulple E values fall within the buffer, then the smallest E is used. Defaults to 0.99, i.e. any value within 99\% of the "best" E can be chosen.
#' @param predtype A character indicating which metric to be used for determining predictive skill. Defaults to "rho". Can be "rho", "mae", or "rmse".
#' @keywords rEDM, embedding dimension
#' @return Value of E corresponding to the best predictive skill.
#' @import stats
#' @export

find_max_theta<-function(smap_out, buffer=0.99, predtype="rho") {

  if(predtype=="rho") {
    kpps<-which(smap_out[,predtype]>=(max(smap_out[,predtype])*buffer))
    best_theta<-smap_out$theta[min(kpps)]
  } else if(predtype%in%c("mae", "rmse")) {
    kpps<-which(smap_out[,predtype]<=(min(smap_out[,predtype])*((1-buffer)+1)))
    best_theta<-smap_out$theta[min(kpps)]
  } else {
    return("error: predtype must be 'rho', 'mae', or 'rmse'")
  }

  return(best_theta)
}
