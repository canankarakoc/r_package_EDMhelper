#' ccm_easy
#'
#' Performs convergent cross mapping using simplex projection across all pairwise comparisons.
#' @param df A data.frame where each column is a time series. Note that additional columns such as time and replicates should not be included in this data.frame. 
#' @param best_E The embedding dimensions to use for time delay embedding.
#' @param lib_segments Segments which are used for predictions.
#' @param pred_segments Segments which are used for predictions  
#' @param ... Additional arguments to be passed to ccm
#' @keywords rEDM, ccm
#' @return A data.frame with forecast statistics
#' @import rEDM
#' @export

ccm_easy <- function (df, best_E, lib_segments, pred_segments=lib_segments, ...){
  
  vars            <- expand.grid(lib_column = colnames(df), target_column=colnames(df))
  vars$best_E     <- unname(best_E[match(vars[,"lib_column"],names(best_E))])
  vars_matrix     <- vars[vars$lib_column!= vars$target_column,]
  vars_matrix_ord <- vars_matrix[order(vars_matrix$lib_column, vars_matrix$target_column),]
  
  
 ccm_output <- do.call(rbind, lapply(seq_len(NROW(vars_matrix_ord)), function(i){
    ccm(df, 
        lib = lib_segments, 
        pred = pred_segments, 
        E=vars_matrix_ord$best_E[i],  
        lib_column = as.character(vars_matrix_ord$lib_column[i]),
        target_column = as.character(vars_matrix_ord$target_column[i]),
        silent=TRUE,...)
    }))

}

