#load data
require(rEDM)
data("paramecium_didinium")
data<-paramecium_didinium
dnames<-colnames(data)[-1]

#get segments (all observations are in one single segment)
sgm<-segment_data(rep(1, nrow(data)))

#make standardized data
df<-normalize_data(data[,dnames])$xnorm

#run simplex on all columns
maxE_test<-round(sqrt(min(sgm[,2]-sgm[,1]))) #maximum E to test is sqrt(time series length)

simplex_out <- lapply(dnames, function(dnames){
  simplex(df[,c(dnames)], E=2:maxE_test, lib=sgm, pred=sgm, tau=1, silent = TRUE)
})

#find best E for all time series
best_E<-find_max_E(simplex_out, buffer=1, predtype = "rho")
names(best_E)<-dnames

ccm_out<-ccm_easy(df=df, best_E=best_E, lib_segments=sgm, pred_segments=sgm)

summaryout<-ccm_summary(ccm_output = ccm_out, predtype = "rho")

sigout<-significance_test(ccm_output = ccm_out, predtype = "rho")
sigout

ccm_easy_plot(summaryout)
