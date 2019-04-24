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
names(simplex_out)<-dnames

plot_simplex(simplex_out)

#find best E for all time series
best_E<-find_max_E(simplex_out, buffer=1, predtype = "rho")
names(best_E)<-dnames

#run CCM
ccm_out<-ccm_easy(df=df, best_E=best_E, lib_segments=sgm, pred_segments=sgm)

#summarize results
summaryout<-ccm_summary(ccm_output = ccm_out, predtype = "rho")

#run significance tests
sigout<-significance_test(ccm_output = ccm_out, predtype = "rho")
sigout

#plot result
ccm_easy_plot(summaryout)

#run time-delay CCM
ccm_del_out<-ccm_delay_easy(df=df, best_E = best_E, lib_segments=sgm, pred_segments=sgm, tp = -5:5)

#summarize results
summaryout_del<-ccm_summary(ccm_output = ccm_del_out, predtype = "rho")

#run significance tests
sigout_del<-significance_test(ccm_output = ccm_del_out, predtype = "rho")
sigout_del

#plot result
ccm_easy_plot(summaryout_del)

#plot prediction skill at the maximum library size as a function of tp
ccm_delay_plot(summaryout_del)

#run s-mapping
smap_out <- lapply(dnames, function(dnames){
  s_map(df[,c(dnames)], E=best_E[dnames], lib=sgm, pred=sgm, tau=1, silent = TRUE)
})
names(smap_out)<-dnames

#find best theta for all time series
best_theta<-find_max_theta(smap_out, buffer=1, predtype = "rho")
names(best_theta)<-dnames

plot_smap(smap_out)



