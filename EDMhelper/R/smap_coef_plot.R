#' smap_coef_plot
#'
#' Plot s-map coefficients.
#' @param index_sp Index number for the desired plot - Calls the corresponding analysis number from smap_coef_out.
#' @param smap_coef_out output of a call "get_smap_coef".
#' @param covar Covariate to plot on x-axis. Defaults to "time". Can be either a column name from the block, or a vector.
#' @param ycol S-map coefficient columns to plot. Defaults to 'all' - i.e. all variables.
#' @param xlab Name to plot on x-axis. Defaults to "".
#' @param ... Additional arguments to be passed to ggplot2.
#' @keywords rEDM, ggplot2
#' @return A ggplot2 object
#' @import rEDM
#' @import ggplot2
#' @import reshape
#' @export

smap_coef_plot <- function (index_sp, smap_coef_out, covar="time", ycol="all", xlab=""){


  colnames(smap_coef_out$smap_out[[index_sp]]$smap_coefficients[[1]])<-c(colnames(smap_coef_out$block_out[[index_sp]])[-c(1:2)], "intercept")


 if(ycol=="all") {
   yvar<-smap_coef_out$smap_out[[index_sp]]$smap_coefficients[[1]]
 } else {
   yvar<-smap_coef_out$smap_out[[index_sp]]$smap_coefficients[[1]][,ycol]
 }
  if(is.null(dim(yvar))) {
    yvar<-as.matrix(yvar)
    colnames(yvar)<-ycol
  }

 if(is.character(covar)) {
   covar<-smap_coef_out$block_out[[index_sp]][,covar]
 }


  indexmat<-t(matrix(nrow=2, rep(which(colnames(smap_coef_out$smap_out[[index_sp]]$smap_coefficients[[1]])%in%colnames(yvar)), each=2)))
  sddat<-sqrt(t(matrix(nrow=ncol(yvar), unlist(lapply(smap_coef_out$smap_out[[index_sp]]$smap_coefficient_covariances[[1]],
                                                      function(x) {
                                                        if(is.null(x)) {
                                                          rep(NA, nrow(indexmat))
                                                        } else {
                                                          x[indexmat]
                                                        }
                                                      })))))


  lower <- suppressMessages(melt(yvar-sddat))
  upper <- suppressMessages(melt(yvar+sddat))

  #Data for plotting
  smap_coef_out_long <- suppressMessages(melt(as.data.frame(yvar)))


  nice_theme<-theme_bw()+
    theme(axis.text= element_text(size=12),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12,face="bold"),
          plot.title = element_text(size=14, face="bold"),
          strip.text = element_text(size=14, face="bold"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())


  value    <- smap_coef_out_long$value
  variable <- smap_coef_out_long$variable

  covar_long<-data.frame(covar=covar, value=value)

  plot_smap_coef <- ggplot(smap_coef_out_long, aes(x=covar_long$covar, y=value, group=variable))+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_ribbon(aes(ymin=lower$value, ymax=upper$value, x=covar_long$covar),alpha=0.1)+
    geom_line(aes(color=variable))+
    #facet_grid(~index_sp)+
    labs(y="S-map coefficients", x=xlab)+
    ggtitle(paste("target: ", smap_coef_out$direction[index_sp,]["target"], sep=""))+
    nice_theme

  return(plot_smap_coef)

}

