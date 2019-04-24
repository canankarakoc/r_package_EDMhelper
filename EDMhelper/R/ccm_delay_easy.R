#' ccm_delaye_easy
#'
#' Performs convergent cross mapping using simplex projection across all pairwise comparisons and a range of prediction lag.
#' @param df A data.frame where each column is a time series. Note that additional columns such as time and replicates should not be included in this data.frame.
#' @param best_E The embedding dimensions to use for time delay embedding.
#' @param lib_segments Segments which are used for predictions.
#' @param pred_segments Segments which are used for predictions.
#' @param tp Range of prediction lag should be tested. Default is -10 to 10.
#' @param ... Additional arguments to be passed to ccm.
#' @keywords rEDM, ccm
#' @return A data.frame with forecast statistics for each prediction lag tested.
#' @import rEDM
#' @export

ccm_delay_easy <- function (df, best_E, lib_segments, pred_segments=lib_segments, tp=-10:10, ...){

  vars            <- expand.grid(lib_column = colnames(df), target_column=colnames(df), tp=tp)
  vars$best_E     <- unname(best_E[match(vars[,"lib_column"],names(best_E))])
  vars_matrix     <- vars[vars$lib_column!= vars$target_column,]
  vars_matrix_ord <- vars_matrix[order(vars_matrix$lib_column, vars_matrix$target_column),]


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


























params   <- expand.grid(lib_column = vars, target_column=vars, tp=-10:10)
params$E <- best_E
params   <- params[params$lib_column != params$target_column,]

ccm_tp_output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i){
  ccm(com_ts,
      lib = segments,
      pred = segments,
      E=params$E[i],
      lib_column = params$lib_column[i],
      target_column = params$target_column[i],
      lib_sizes = NROW(data_by_rep[[1]]),
      tp=params$tp[i],
      random_libs = T,
      num_samples = 100, #10000 saved
      silent=T)
}))

#saveRDS(ccm_tp_output, "CCM_delayed_predation.R")
ccm_tp_output <- readRDS("CCM_delayed_predation.R")

ccm_tp_output$direction <- paste(ccm_tp_output$lib_column, "xmap to", ccm_tp_output$target_column)

#CCM_Means
mylist <- split(ccm_tp_output, list(ccm_tp_output$direction, ccm_tp_output$tp))
means <- do.call(rbind, lapply(mylist, function(x) ccm_means(x)))
means$direction <- paste(means$lib_column, "xmap to", means$target_column)
