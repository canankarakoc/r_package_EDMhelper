#' get_smap_coef
#'
#' Calculate s-mapping coefficients based on results from ccm-test
#' @param df A data.frame where each column is a time series. Note that additional columns such as time and replicates should not be included in this data.frame.
#' @param lib_segments Segments which are used for predictions.
#' @param sigout results of ccm significance test, from significance_test() or ccm_summary() functions
#' @param best_E vector of best E values from find_max_E()
#' @param best_theta vector of best theta values from find_max_theta()
#' @param selfref Include separate analysis of self-limitation? Defaults to FALSE.
#' @param ... Additional arguments to be passed to block_lnlp - e.g. if different pred segments are desired
#' @keywords rEDM, s-mapping
#' @return List including the direction of tets, list of s-mapping results, and list of blocks of data corresponding to the tests.
#' @export

get_smap_coef<-function(df, lib_segments, sigout, best_E, best_theta, selfref=FALSE, ...) {
  direction<-t(matrix(nrow=2, data=unlist(strsplit(as.character(sigout$direction), " causes ", fixed=T))))
  colnames(direction)<-c("target", "lib")

  #which lags to use for reconstructions?
  tpuse<-sigout$tp
  tpuse[tpuse>0]<-0 #future states cannot influence past states

  #add self-limitation
  splst<-sort(unique(c(direction)))
  if(selfref) {
    direction<-rbind(direction, t(matrix(nrow=2, data=rep(splst, each=2))))
    tpuse<-c(tpuse, rep(0, length(splst)))
  }

  smap_out<-list()
  block_out<-list()

  for(i in 1:nrow(direction)) {
    target<-direction[i,"target"]
    lib<-direction[i,"lib"]

    #build block for lib, with one lag replaced by target
    if(lib!=target) {
      block<-make_block(df[,lib], lib = lib_segments, max_lag = pmax(1, best_E[which(names(best_E)==lib)]-1),restrict_to_lib = TRUE)
    } else {
      block<-make_block(df[,lib], lib = lib_segments, max_lag = pmax(1, best_E[which(names(best_E)==lib)]),restrict_to_lib = TRUE)
    }

    if(target==lib) {
      block_target<-NULL
    } else {
      block_target<-make_block(df[,target], lib = lib_segments, max_lag = pmax(1, -tpuse[i]),restrict_to_lib = TRUE)
      block_target<-block_target[,ncol(block_target)]
    }

    if(!is.null(block_target)) {
      block<-cbind(block[,1:2], block_target, block[,-c(1:2)])
    }

    #rename block columns
    if(ncol(block)>3) {
      colnames(block)<-c("time",
                       paste(lib, paste("_t", 1, sep=""), sep=""),
                       paste(target, "_t", abs(tpuse[i]), sep=""),
                       paste(lib, paste("_t", (1:(ncol(block)-2))[-1], sep=""), sep=""))
    } else {
      colnames(block)<-c("time",
                         paste(lib, paste("_t", 1, sep=""), sep=""),
                         paste(target, "_t", abs(tpuse[i]), sep=""))
    }

    smap_out_tmp<-block_lnlp(block=block, method = "s-map",
                             tp = 0, target_column = 1, theta = unlist(best_theta[which(names(best_theta)==lib)]), first_column_time = TRUE, columns = (1:(ncol(block)-1))[-1],
                             stats_only = FALSE, silent = TRUE, save_smap_coefficients = TRUE, ...)

    #note - error in code below - block already removes elements that span a gap
    #smap_out_tmp<-block_lnlp(block=block, lib = lib_segments, method = "s-map",
    #           tp = 0, target_column = 1, theta = unlist(best_theta[which(names(best_theta)==lib)]), first_column_time = TRUE, columns = (1:(ncol(block)-1))[-1],
    #           stats_only = FALSE, silent = TRUE, save_smap_coefficients = TRUE, ...)

    smap_out[[i]]<-smap_out_tmp
    block_out[[i]]<-block

    colnames(smap_out[[i]]$smap_coefficients[[1]])<-c(colnames(block_out[[i]])[-c(1:2)], "intercept")
  }

  return(list(direction=direction, smap_out=smap_out, block_out=block_out))
}
