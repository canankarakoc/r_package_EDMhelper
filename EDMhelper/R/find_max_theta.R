#' find_max_theta
#'
#' Finds the best nonlinearity parameter theta corresponding to the best predictive skill for a simplex call.
#' @param smap_out output of a call from "smap"
#' @param buffer A number from 0 to 1, indicating the allowed buffer distance between the "best" theta, and the theta that is used. If multiple values fall within the buffer, then the smallest theta is used. Defaults to 0.99, i.e. any value within 99\% of the "best" theta can be chosen.
#' @param predtype A character indicating which metric to be used for determining predictive skill. Defaults to "rho". Can be "rho", "mae", or "rmse".
#' @keywords rEDM, nonlinearity parameter
#' @return Value of theta corresponding to the best predictive skill.
#' @import stats
#' @export

find_max_theta<-function(smap_out, buffer=0.99, predtype="rho") {

  if(sum(grep("const_p_val", names(smap_out)))==0) {
    best_theta<-lapply(smap_out, function(x) find_max_theta(x, buffer=buffer, predtype=predtype))
  } else {
    if(predtype=="rho") {
      kpps<-which(smap_out[,predtype]>=(max(smap_out[,predtype])*buffer))
      best_theta<-smap_out$theta[min(kpps)]
    } else if(predtype%in%c("mae", "rmse")) {
      kpps<-which(smap_out[,predtype]<=(min(smap_out[,predtype])*((1-buffer)+1)))
      best_theta<-smap_out$theta[min(kpps)]
    } else {
      return("error: predtype must be 'rho', 'mae', or 'rmse'")
    }
  }

  return(best_theta)
}
