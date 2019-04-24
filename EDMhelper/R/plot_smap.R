#' plot_smap
#'
#' Plot prediction skill vs. theta for the results of an s-map
#' @param smap_out Output from s-map function (or a list of outputs)
#' @param predtype Type of prediction metric to use. Can be "rho", "rmse", or "mae".
#' @param pname Optional name for the plot
#' @keywords rEDM, s-map
#' @return A ggplot object
#' @import ggplot2
#' @export

plot_smap<-function(smap_out, predtype="rho", pname=NA) {

  if(sum(grep("const_p_val", names(smap_out)))==0) {
    dm<-dim(smap_out[[1]])
    lng<-length(smap_out)
    if(!is.na(pname[1])) {
      nms<-pname
    } else {
      nms<-names(smap_out)
    }
    smap_out<-do.call("rbind", smap_out)
    rownames(smap_out)<-NULL
    smap_out$variable<-rep(nms, each=dm[1])

  } else {
    smap_out$variable<-pname
  }
  variable<-smap_out$variable

  theta<-smap_out$theta


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

  predtype<-smap_out[,predtype]

  plot_ccm <- ggplot(smap_out, aes(x=theta, y=predtype))+
    geom_line(color="darkred")+
    facet_wrap(~variable)+
    labs(y="Prediction skill", x="Nonlinearity, Theta")+
    nice_theme

  return(plot_ccm)
}
