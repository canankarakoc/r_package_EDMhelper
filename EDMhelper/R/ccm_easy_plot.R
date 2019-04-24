#' ccm_easy_plot
#'
#' Plots convergent cross mapping results.
#' @param ccm_summary output of a call "ccm_easy".
#' @param ... Additional arguments to be passed to ggplot2.
#' @keywords rEDM, ggplot2
#' @return A ggplot2 object
#' @import rEDM, ggplot2
#' @export

ccm_easy_plot <- function (cmm_summary){


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



 plot_ccm <- ggplot(ccm_summary, aes(x=lib_size, y=Q...))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_ribbon(aes(ymin=Q..., ymax=Q..., x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  facet_wrap(~direction)+
  labs(y="Prediction skill", x="Library size")+
  nice_theme

  return(plot_ccm)

}

