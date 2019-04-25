#' ccm_delay_plot
#'
#' Plot prediction skills as a function of prediction lag (tp).
#' @param summaryout output of a call "ccm_summary".
#' @param ... Additional arguments to be passed to ggplot2.
#' @keywords rEDM, ggplot2
#' @return A ggplot2 object.
#' @import rEDM
#' @import ggplot2
#' @export

ccm_delay_plot <- function (summaryout){

  nice_theme <-theme_bw()+
    theme(axis.text= element_text(size=12),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12,face="bold"),
          plot.title = element_text(size=14, face="bold"),
          strip.text = element_text(size=14, face="bold"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  summaryout_maxminlibsize <-  summaryout[which(summaryout$lib_size==max(summaryout$lib_size)
                                          | summaryout$lib_size==min(summaryout$lib_size)),]


  tp         <- summaryout_maxminlibsize$tp
  lib_size   <- summaryout_maxminlibsize$lib_size
  Q50.       <- summaryout_maxminlibsize$Q50.
  Q15.86553. <- summaryout_maxminlibsize$Q15.86553.
  Q84.13447. <- summaryout_maxminlibsize$Q84.13447.


  plot_ccm_delay <- ggplot(summaryout_maxminlibsize, aes(x=tp, y=Q50., group=as.factor(lib_size)))+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_vline(xintercept=0, linetype="dashed")+
    geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447.),alpha=0.1)+
    geom_line(aes(color=as.factor(lib_size)))+
    facet_grid(~direction)+
    labs(y="Prediction skill", x="tp")+
    scale_x_continuous(breaks = tp)+
    scale_color_manual(name="Lib. size", values=c("darkblue", "darkred"))+
    nice_theme

  return(plot_ccm_delay)

}

