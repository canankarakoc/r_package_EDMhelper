get_smap_coef<-function(sigout, best_E, best_theta) {
  direction<-t(matrix(nrow=2, data=unlist(strsplit(as.character(sigout$direction), " causes ", fixed=T))))
  colnames(direction)<-c("target", "lib")
  
  
  
}
