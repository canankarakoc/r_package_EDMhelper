#load data
require(rEDM)
data("sardine_anchovy_sst")
dnames<-c("anchovy", "sardine", "sio_sst", "np_sst")

#get segments (all observations are in one single segment)
sgm<-segment_data(rep(1, nrow(sardine_anchovy_sst)))

#make standardized data
df<-normalize_data(sardine_anchovy_sst[,dnames])$xnorm

#run simplex on all columns
maxE_test<-round(sqrt(min(sgm[,2]-sgm[,1]))) #maximum E to test is sqrt(time series length)
simplex_out <- lapply(dnames, function(dnames){
  simplex(df[,c(dnames)], E=2:maxE_test, lib=sgm, pred=sgm, tau=1, silent = TRUE)
})

#find best E for all time series
best_E <- unlist(lapply(simplex_out, function(x) find_max_E(x, buffer = 1, predtype = "rho")))
names(best_E)<-dnames

ccm_out<-ccm_easy(df=df, best_E=best_E, lib_segments=sgm, pred_segments=sgm)

summaryout<-ccm_summary(ccm_out)
