#' ccm_delay_plot
#'
#' Plots prediction skills as a function of prediction lag tp.
#' @param summaryout output of a call "ccm_summary".
#' @param ... Additional arguments to be passed to ggplot2.
#' @keywords rEDM, ggplot2
#' @return A ggplot2 object.
#' @import rEDM
#' @import ggplot2
#' @export

ccm_delay_plot <- function (summaryout){

  nice_theme<-theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title= element_text(size=12,face="bold"),
          plot.title = element_text(size=14, face="bold"),
          strip.text =element_text(size=12),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  summaryout_maxlibsize <-  summaryout[which(summaryout$lib_size==max(summaryout$lib_size)),]

  plot_ccm_delay <- ggplot(summaryout_maxlibsize, aes(x=tp, y=Q50.))+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447.),alpha=0.1)+
    geom_line(color="darkred")+
    facet_grid(~direction)+
    labs(y="Prediction skill", x="tp")+
    nice_theme

  return(plot_ccm_delay)

}

