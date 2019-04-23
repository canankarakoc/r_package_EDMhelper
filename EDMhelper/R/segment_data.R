#' segment_data
#'
#' Splits data by replicates, create segments for predictions
#' @param replicates A data column contains the replicates to be splitted
#' @keywords split, replicates
#' @return A matrix including the beginning and the end of the segments
#' @export
#' @examples
#' fakedata = cbind(rep(1:3, each=10), rnorm(30))
#' segment_data(fakedata[,1])


segment_data<-function(replicates) {
  
  segments_end <-  cumsum(table(replicates))
  segments_begin <- c(1, segments_end[-length(segments_end)]+1)
  segmentout = cbind(segments_begin,segments_end)
  rownames(segmentout) = NULL
  return(segmentout)
}


