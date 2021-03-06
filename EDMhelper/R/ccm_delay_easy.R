#' ccm_delay_easy
#'
#' Performs convergent cross mapping using simplex projection across all pairwise comparisons and a range of prediction lag.
#' @param df A data.frame where each column is a time series. Note that additional columns such as time and replicates should not be included in this data.frame.
#' @param best_E The embedding dimensions to use for time delay embedding.
#' @param lib_segments Segments which are used for predictions.
#' @param pred_segments Segments which are used for predictions.
#' @param tp Range of prediction lag should be tested. Default is -10 to 10.
#' @param vars_matrix_ord An optional pre-made data.frame of combinations of variables to test. Must include columns 'lib_column', 'target_column', and 'best_E'.
#' @param ... Additional arguments to be passed to ccm.
#' @keywords rEDM, ccm
#' @return A data.frame with forecast statistics for each prediction lag tested.
#' @import rEDM
#' @export

ccm_delay_easy <- function (df, best_E, lib_segments, pred_segments=lib_segments, tp=-3:3, vars_matrix_ord=NULL, ...){

  if(is.null(vars_matrix_ord)) {
    vars            <- expand.grid(lib_column = colnames(df), target_column=colnames(df), tp=tp)
    vars$best_E     <- unname(best_E[match(vars[,"lib_column"],names(best_E))])
    vars_matrix     <- vars[vars$lib_column!= vars$target_column,]
    vars_matrix_ord <- vars_matrix[order(vars_matrix$lib_column, vars_matrix$target_column),]
  } else {
    vars_matrix_ord<-data.frame(vars_matrix_ord[rep(1:nrow(vars_matrix_ord), each=length(tp)),], tp=rep(tp, nrow(vars_matrix_ord)))
    vars_matrix_ord<-vars_matrix_ord[,c("lib_column", "target_column", "tp", "best_E")]
  }

  ccm_delay_output <- do.call(rbind, lapply(seq_len(NROW(vars_matrix_ord)), function(i){
    ccm(df,
          lib = lib_segments,
          pred = pred_segments,
          E=vars_matrix_ord$best_E[i],
          lib_column = as.character(vars_matrix_ord$lib_column[i]),
          target_column = as.character(vars_matrix_ord$target_column[i]),
          tp=vars_matrix_ord$tp[i],
          silent=TRUE,...)
  }))

  return(ccm_delay_output)
}
