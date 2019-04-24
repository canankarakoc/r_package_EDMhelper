#' plot_simplex
#'
#' Plot prediction skill vs. theta for the results of the simplex function
#' @param simplex_out Output from simplex function (or a list of outputs)
#' @param predtype Type of prediction metric to use. Can be "rho", "rmse", or "mae".
#' @param pname Optional name for the plot
#' @keywords rEDM, s-map
#' @return A ggplot object
#' @import ggplot2
#' @export

plot_simplex<-function(simplex_out, predtype="rho", pname=NA) {
  
  if(sum(grep("const_p_val", names(simplex_out)))==0) {
    dm<-dim(simplex_out[[1]])
    lng<-length(simplex_out)
    if(!is.na(pname[1])) {
      nms<-pname
    } else {
      nms<-names(simplex_out)
    }
    simplex_out<-do.call("rbind", simplex_out)
    rownames(simplex_out)<-NULL
    simplex_out$variable<-rep(nms, each=dm[1])
      
  } else {
    simplex_out$variable<-pname
  }
  variable<-simplex_out$variable
  
  E<-simplex_out$E
  
  
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
  
  predtype<-simplex_out[,predtype]
  
  plot_ccm <- ggplot(simplex_out, aes(x=E, y=predtype))+
    geom_line(color="darkred")+
    facet_wrap(~variable)+
    labs(y="Prediction skill", x="Embedding dimension, E")+
    nice_theme
    
  return(plot_ccm)
}
