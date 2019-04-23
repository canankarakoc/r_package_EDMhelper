#' normalize_data
#'
#' Normalizes data, either for data with or without true zeros.
#' @param x a numeric vector to be normalized
#' @param truezero A logical telling whether the data includes true zeros - if TRUE then standardizes based only on the standard deviaiton of x, else if FALSE standardizes x based on its mean and standard deviation. Defaults to FALSE.
#' @keywords normalization
#' @return A list, including the standardized values of x (xnorm), and the mean (mux) and the standard deviation (sdx) used for the transformation.
#' @export

normalize_data<-function(x, truezero=FALSE) {
  
  sdx<-sd(x,na.rm=TRUE)
  if(truezero) {
    mux<-0
  } else {
    mux<-mean(x, na.rm=T)
  }
  xnorm<-(x-mux)/sdx
  
  return(list(xnorm=xnorm, mux=mux, sdx=sdx))
}

#' inormalize_data
#'
#' Inverse transformation of normalize_data function.
#' @param normout Output from normalize_data function.
#' @param x An optional vector of values to be untranformed based on the mean and standard deviation stored in normalize_data. Defaults to NA, in which case the values in normout are used.
#' @keywords normalization
#' @return Untransformed values
#' @export

inormalize_data<-function(normout, x=NA) {
  
  if(is.na(x)) {
    x<-normout$xnorm
  }
  x_untr<-x*normout$sdx+normout$mux
  
  return(x_untr)
}
