#' stability_analysis
#'
#' Calculate local stability and interaction indices. Follows methods of Ushio et al. 2018.
#' @param smap_coef_out Output of a call "get_smap_coef".
#' @param doplot Logical argument indicating if indices should be plotted or not. Default is TRUE.
#' @keywords rEDM, s-mapping
#' @return List including indices, and plots if requested.
#' @source Usho, M., et al. (2018). Fluctuating interaction network and time-varying stability of a natural fish community. Nature 554:360–363.
#' @import ggplot2
#' @export
#'

stability_analysis <- function (smap_coef_out, doplot=TRUE){

  newmat_tot = NULL
  for (i in 1:nrow(smap_coef_out$direction)){


    sppos <- i

    tmp = (smap_coef_out$smap_out[[sppos]]$smap_coefficients[[1]])

    #self      <-  smap_coef_out$smap_out[[sppos]]$smap_coefficients[[1]][,1, drop=FALSE]
    if(smap_coef_out$direction[i,"target"]!=smap_coef_out$direction[i,"lib"]) {
      variables <-  tmp[,grep(smap_coef_out$direction[i,"target"], colnames(tmp)), drop=FALSE]
    } else {
      variables <-  tmp[,grep("0", colnames(tmp)), drop=FALSE]
    }

    newmatrix <- cbind(variables,
                       target=rep(unname(smap_coef_out$direction[sppos,'target'])),
                       lib=rep(unname(smap_coef_out$direction[sppos,'lib'])),
                       time = smap_coef_out$block_out[[sppos]][,'time'])
    colnames(newmatrix)[1] <- 'value'

    newmat_tot = rbind(newmat_tot, newmatrix)

  }

 newmat_tot_sorted         <-  newmat_tot[order( newmat_tot$time),]

 newmat_tot_sorted_splited <- split(newmat_tot_sorted, newmat_tot_sorted$time)
 totlst = sort(unique(c(as.character(newmat_tot$target), as.character(newmat_tot$lib))))

 spmatout =  lapply(newmat_tot_sorted_splited, function(x) {
    psmat = cbind(match(x$target, totlst), match(x$lib, totlst))
    matout = matrix(nrow = length(totlst), ncol=length(totlst), data=0)
    matout[psmat] = x$value
    matout
  })

 #Local stability
 local_stability <- sapply(spmatout, function(x){
   if(all(is.finite(x))) {
     max(Re(eigen(x)$values))
   } else {
     NA
   }
   })

 #Interaction indexes
 #Mean interactions index
 mean_interactions <- sapply(spmatout, function(x){mean(abs(x[row(x)!=col(x)]))})

 #Weak interactions index
 time_median <- sapply(spmatout, function(x){median(abs(x[row(x)!=col(x)]))})
 time_max <- sapply(spmatout, function(x){max(abs(x[row(x)!=col(x)]))})
 weak_interactions <- as.vector(time_median)/as.vector(time_max)

 time = smap_coef_out$block_out[[1]][,'time']

 all_indices <- data.frame(time, local_stability, mean_interactions,
                           weak_interactions)



 if (doplot){

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


  all_indices_melted <- melt(all_indices, id="time")
  variable           <- all_indices_melted$variable
  value              <- all_indices_melted$value



   plot_indices <- ggplot(all_indices_melted, aes(x=time, y=value))+
     #geom_hline(yintercept=0, linetype="dashed")+
     geom_line(color="darkred")+
     facet_wrap(~variable, scales="free")+
     labs(y="Index", x="Time")+
     nice_theme


 }else{
   plot_indices <- NULL
 }

 return(list(all_indices=all_indices, plot_indices=plot_indices))


}


