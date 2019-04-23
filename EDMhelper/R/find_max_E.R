#' find_max_E
#'
#' Finds the embedding dimension corresponding to the best predictive skill for a simplex call.
#' @param simplex_out output of a call from "simplex"
#' @param buffer A number from 0 to 1, indicating the allowed buffer distance between the "best" E, and the E that is used. If mulple E values fall within the buffer, then the smallest E is used. Defaults to 0.99 (i.e. any value within 99% of the "best" E can be chosen)
#' @param type A character indicating which metric to be used for determining predictive skill. Defaults to "rho". Can be "rho", "mae", or "rmse".
#' @keywords rEDM, embedding dimension
#' @return Value of E corresponding to the best predictive skill.
#' @export
#' @import

find_max_E<-function(simplex_out, buffer=0.99, type="rho") {

  if(type=="rho") {
    kpps<-which(simplex_out[,type]>=(max(simplex_out[,type])*buffer))
    best_E<-simplex_out$E[min(kpps)]
  } else if(type%in%c("rmse", "rmse")) {
    kpps<-which(simplex_out[,type]<=(min(simplex_out[,type])*((1-buffer)+1)))
    best_E<-simplex_out$E[min(kpps)]
  } else {
    return("error: type must be 'rho', 'mae', or 'rmse'")
  }

  return(best_E)
}
