#' ccm_easy_plot
#'
#' Plots convergent cross mapping results.
#' @param summaryout output of a call "ccm_easy".
#' @param ... Additional arguments to be passed to ggplot2.
#' @keywords rEDM, ggplot2
#' @return A ggplot2 object
#' @import rEDM
#' @import ggplot2
#' @export

ccm_easy_plot <- function (summaryout){

  lib_size<-summaryout$lib_size

  nice_theme<-theme_bw()+
    theme(axis.text= element_text(size=12),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12,face="bold"),
          plot.title = element_text(size=12, face="bold"),
          strip.text = element_text(size=12, face="bold"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  Q50.<-summaryout$Q50.
  Q15.86553.<-summaryout$Q15.86553.
  Q84.13447.<-summaryout$Q84.13447.

  plot_ccm <- ggplot(summaryout, aes(x=lib_size, y=Q50.))+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447., x=lib_size),alpha=0.1)+
    geom_line(color="darkred")+
    facet_grid(tp~direction)+
    labs(y="Prediction skill", x="Library size")+
    nice_theme

  return(plot_ccm)

}

