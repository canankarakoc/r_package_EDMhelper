#' normalize_data
#'
#' Normalizes data, either for data with or without true zeros.
#' @param x a numeric vector to be normalized, or a matrix where each column is a time series
#' @param truezero A logical telling whether the data includes true zeros - if TRUE then standardizes based only on the standard deviaiton of x, else if FALSE standardizes x based on its mean and standard deviation. Defaults to FALSE.
#' @keywords normalization
#' @return A list, including the standardized values of x (xnorm), and the mean (mux) and the standard deviation (sdx) used for the transformation.
#' @export

normalize_data<-function(x, truezero=FALSE) {

  if(length(dim(x))==2) {
    nps<-ncol(x)
  } else {
    nps<-1
    x<-as.matrix(x)
    if(ncol(x)>1) {
      x<-t(x)
    }
  }

  xnorm<-matrix(nrow=nrow(x), ncol=nps)
  mux<-numeric(length(nps))
  sdx<-numeric(length(nps))

  for(i in 1:nps) {
    sdx[i]<-sd(x[,i],na.rm=TRUE)
    if(truezero) {
      mux[i]<-0
    } else {
      mux[i]<-mean(x[,i], na.rm=T)
    }
    xnorm[,i]<-(x[,i]-mux[i])/sdx[i]
  }
  colnames(xnorm)<-colnames(x)

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

  if(length(dim(x))==2) {
    nps<-ncol(x)
  } else {
    nps<-1
    x<-as.matrix(x)
    if(ncol(x)>1) {
      x<-t(x)
    }
  }

  x_untr<-matrix(nrow=nrow(x), ncol=nps)
  for(i in 1:nps) {
    x_untr[,i]<-x[,i]*normout$sdx[i]+normout$mux[i]
  }
  colnames(x_untr)<-colnames(x)

  return(x_untr)
}
