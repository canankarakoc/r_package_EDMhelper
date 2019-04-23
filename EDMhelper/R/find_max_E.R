#' find_max_E
#'
#' Finds the embedding dimension corresponding to the best predictive skill for a simplex call.
#' @param simplex_out output of a call from "simplex"
#' @param df a data frame that includes columns "E" (embedding dimension) and "rho" (prediction skill)
#' @keywords rEDM, embedding dimension
#' @return Value of E corresponding to the best predictive skill.
#' @export
#' @import rEDM

#Here I'd write a function max E finds minimum E which maximizes 10 percent. 
#I am not sure if it is stupid, I rounded rhos

find_max_E<-function(simplex_out, df) {
  max_E <- sapply(simplex_out, function(df){
    df$E[which.max(round(df$rho,1))]
  })
  max_E
}
