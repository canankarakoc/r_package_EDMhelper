#Title: Increase in causal forcings in a predator-prey system
#Authors: CK, ATC
#Date: 7/11/19 - 27/2/19

library(rEDM)
library(reshape)
library(ggplot2)
#library(cowplot)
library(plyr)
#library(zoo)
#library(dplyr)# cause problems with plyr
library(tidyr)
#library(geomnet)
#library(vegan)


mytheme<-theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        legend.text = element_text(size=12), 
        legend.title= element_text(size=12,face="bold"), 
        plot.title = element_text(size=14, face="bold"), 
        strip.text =element_text(size=12))

#Time series of prey species with/without predator
#Data comes from microcosms which were maintained for 34 days by 10% daily transfer 
#of the cultures into fresh media. Abundances of all the prey species were equal 
#at the start. Counting started on day three. All other culture conditions were constant. 

#There is autocorrelation in time. I took the first difference before ccm but then 
#rho decreases dramatically.

setwd("~/Everyfthing/a-CCM-Granger/EDM")
source("Summary_statistics_functions.R")
predation.data   <- read.table("p_combined.csv", header=T,sep=";",dec=",", stringsAsFactors=F)
competition.data <- read.table("wp_combined.csv", header=T,sep=";",dec=",", stringsAsFactors=F)

str(predation.data)
str(competition.data)

#Fill missing values with linear intrapolation 
predation <- as.data.frame(cbind(day=predation.data$day, 
                                 ddply(predation.data[,2:7], 
                                       .(replicate=as.factor(predation.data$replicate)), 
                                       colwise(na.approx))))

competition <- as.data.frame(cbind(day=competition.data$day, 
                                   ddply(competition.data[,c(2:5)], 
                                         .(replicate=as.factor(competition.data$replicate)), 
                                         colwise(na.approx))))


#Plot the time series 
predation_ts <-melt(predation, id=c("day","replicate"))
competition_ts <-melt(competition, id=c("day","replicate"))

dfc_pred  <- summarySE(predation_ts, measurevar="value", groupvars=c("day","variable"),na.rm=T)
dfc_comp  <- summarySE(competition_ts, measurevar="value", groupvars=c("day","variable"),na.rm=T)

#With predator (competition and predation)
pred.mean<- ggplot(dfc_pred, aes(x=day,y=log10(value)))+
  #geom_vline(xintercept=3, color="#FFCC00", size=1.5, alpha=0.6)+
  #geom_segment(mapping=aes(x=3, y=7.2, xend=5, yend=7.2), arrow=arrow(), 
  #             size=1, color="#FFCC00", alpha=.6)+
  geom_errorbar(aes(ymin=log10((value-sd)),ymax=log10((value+sd))), 
                color="grey50", width=0.2, size=0.5)+ 
  geom_line(size=1, aes(color=variable))+
  geom_point(size=2, aes(fill=variable, shape=variable))+
  ylab("log(Abundance)")+ 
  xlab("Time (days)")+
  #ggtitle("Predation")+
  scale_shape_manual(values=c(21,21,21,21,21,24), name = "Species")+
  scale_fill_manual(values=c("grey50", "grey50","midnightblue", "skyblue3", 
                             "grey50", "grey25"), name = "Species")+
  scale_color_manual(values=c("grey50", "grey50","midnightblue", "skyblue3", 
                              "grey50", "grey25"), name = "Species")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
  
pred.mean

#Without predator (only competition)
comp.mean<- ggplot(dfc_comp, aes(x=day,y=log10(value)))+
  #geom_vline(xintercept=3, color="#FFCC00", size=1.5, alpha=0.6)+
  #geom_segment(mapping=aes(x=3, y=9.2, xend=5, yend=9.2), arrow=arrow(), 
  #             size=1, color="#FFCC00", alpha=.6)+
  geom_errorbar(aes(ymin=log10((value-sd)),ymax=log10((value+sd))), 
                color="grey50", width=0.2, size=0.5)+ 
  geom_line(size=1, aes(color=variable))+
  geom_point(size=2, aes(fill=variable, shape=variable))+
  ylab("log(Abundance)")+ 
  xlab("Time (days)")+
  #ggtitle("Competition")+
  scale_shape_manual(values=c(21,21,21,21,21,24), name = "Species")+
  scale_fill_manual(values=c("darkred", "lightpink3",  "midnightblue",
                             "darkgreen"), name = "Species")+
  scale_color_manual(values=c("darkred", "lightpink3", "midnightblue", "darkgreen"), name = "Species")+
  #coord_polar()+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
  
comp.mean

plot_grid(pred.mean, comp.mean, ncol=1,
          align = 'h', rel_widths = c(1, 1))

#T is the predator, others are the prey species, m indicates the different phenotype 
#of a species. Am forms biofilm in the liquid gas interface, Km forms huge clumps which 
#are resistant to predation. 

#*****************************************
#Analysis for the microcosms with predator
#*****************************************
predation_norm <- as.data.frame(apply(predation[,3:8], 2, function(x) (x - min(x))/(max(x)-min(x))))
competition_norm <- as.data.frame(apply(competition[,c(3,4,6)], 2, function(x) (x - min(x))/(max(x)-min(x))))

predation   <- cbind(predation[,1:2], predation_norm)
competition <- cbind(competition[,1:2], competition_norm)

#Separate time column from variables 
#Split data by replicates
#Use all the segments for prediction 

vars <- c("A","Am","K","Km","S","T")
com_ts <- predation[,vars]
data_by_rep <- split(com_ts, predation$replicate)
segments_end <- cumsum(sapply(data_by_rep, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)]+1)
segments <- cbind(segments_begin,segments_end)

#Predictability & nonlinearity
#Find best E
simplex_out <- lapply(vars, function(var){
  simplex(predation[,c("day", var)], E=1:6, lib=segments, pred=segments, tau=1)
})

names(simplex_out) <- vars

par(mfrow=c(3,2))
for(var in names(simplex_out)){
  plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, type="l", 
       xlab="Embedding Dimension (E)", ylab="Forecast Skill (rho)", main = var)
}

#Here I'd write a function max E finds minimum E which maximizes 10 percent. 
#I am not sure if it is stupid, I rounded rhos 
max_E <- sapply(simplex_out, function(df){
  df$E[which.max(round(df$rho,1))]
})
max_E

#Choose the smaller E if the difference in rho is not bigger than ~10%
best_E <- as.integer(c(5,4,2,4,4,6))
names(best_E) <- vars

##Nonlinearity
smap_out <- lapply(vars, function(var){
  s_map(predation[,c("day", var)], E=best_E[var], lib= segments)
})
names(smap_out) <- names(simplex_out)

par(mfrow=c(3,2))
for(var in names(smap_out)){
  plot(smap_out[[var]]$theta, smap_out[[var]]$rho, type="l", 
       xlab="Nonlinearity (theta)", ylab="Forcast Skill (rho)", main=var)
}

#*********
###CCM###
#********

#CCM for tp=0 which is default
params2 <- expand.grid(lib_column = vars, target_column=vars)
params2$E <- best_E
params2 <- params2[params2$lib_column!= params2$target_column,]
params2.ord <-params2[order(params2$lib_column, params2$target_column),]

ccm_output <- do.call(rbind, lapply(seq_len(NROW(params2.ord)), function(i){
  ccm(com_ts, 
      lib = segments, 
      pred = segments, 
      E=params2.ord$E[i],  
      lib_column = params2.ord$lib_column[i],
      target_column = params2.ord$target_column[i],
      random_libs = T,
      num_samples = 100, #10000 saved
      silent=T)
}))

#saveRDS(ccm_output, "CCM_tp0_predation.R")
ccm_output <- readRDS("CCM_tp0_predation.R")

ccm_output$direction <- paste(ccm_output$lib_column, "xmap to", ccm_output$target_column)

ccm_ag <- with(ccm_output, aggregate(cbind(rho=rho),
                                     list(direction=direction,
                                          lib_column=lib_column,
                                          lib_size=lib_size),
                                     function(x) quantile(x, c(0.025,
                                                               pnorm(-1,0,1), 0.5, 
                                                               pnorm(1,0,1), 0.975), 
                                                          na.rm=T)))
#Create a pairs column for plotting 
#Probably so dumm 
ccm_ag$pairs <-mapvalues(ccm_ag$direction, from = c("A xmap to Am", "Am xmap to A", 
                                     "A xmap to K", "K xmap to A",
                                     "A xmap to Km", "Km xmap to A",
                                     "A xmap to S", "S xmap to A",
                                     "A xmap to T", "T xmap to A",
                                     "Am xmap to K", "K xmap to Am",
                                     "Am xmap to Km", "Km xmap to Am",
                                     "Am xmap to S", "S xmap to Am", 
                                     "Am xmap to T", "T xmap to Am",
                                     "K xmap to Km", "Km xmap to K", 
                                     "K xmap to S", "S xmap to K",
                                     "K xmap to T", "T xmap to K", 
                                     "Km xmap to S", "S xmap to Km", 
                                     "Km xmap to T", "T xmap to Km",
                                     "S xmap to T", "T xmap to S"), 
          to=c("A-Am", "A-Am", "A-K", "A-K", "A-Km", "A-Km", "A-S", "A-S", "A-T", "A-T", 
               "Am-K", "Am-K", "Am-Km", "Am-Km", "Am-S", "Am-S", "Am-T", "Am-T", 
               "K-Km","K-Km", "K-S", "K-S", "K-T", "K-T",  "Km-S", "Km-S", "Km-T", 
               "Km-T", "S-T", "S-T"))
                                   

#Significance criteria according to Clark et al. Ecology paper 
func1 <- function(data){1-mean((data[data$lib_size==min(data$lib_size),]$rho<data[data$lib_size==max(data$lib_size),]$rho)&
                                 (data[data$lib_size==max(data$lib_size),]$rho>0), na.rm=T)
}

signif <- ddply(ccm_output,"direction",func1)

#Significance criteria, less conservative
#If terminal rho is bigger than 0 and the difference between terminal rho and initial rho is
#bigger than 0 (convergence)
#In Nature fish paper it'S 0.1

signif$V2 <- ifelse(ccm_ag[ccm_ag$lib_size==max(ccm_ag$lib_size),]$rho[,1]>0 &
                      (ccm_ag[ccm_ag$lib_size==max(ccm_ag$lib_size),]$rho[,1]-
                         ccm_ag[ccm_ag$lib_size==min(ccm_ag$lib_size),]$rho[,1])>0, "1","0")

signif

#Plot 
ccm_ag$guide <- rep(0, length=nrow(ccm_ag))
ccm_fig <- ggplot(ccm_ag, aes(x=lib_size, y=rho[,3], group=lib_column))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  facet_wrap(~pairs)+
  scale_color_manual(name= "Library\ncolumn", values = c("darkred", "lightpink3", "midnightblue",   
                                "skyblue3", "darkgreen", "grey25"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig

#Example plots for presentation 
#A-K
ccm_ag_ak <- subset(ccm_ag, ccm_ag$pairs=="A-K")
ccm_ag_ak$guide <- rep(0, length=nrow(ccm_ag_ak))
ccm_fig_ak <- ggplot(ccm_ag_ak, aes(x=lib_size, y=rho[,3], group=lib_column))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  scale_color_manual(name="Forcing", labels=c("K causes A", "A causes K"), values = c("darkred", "midnightblue"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  theme(legend.position = c(0.8, 0.12), legend.background = element_rect(fill = "transparent"))
ccm_fig_ak

#A-K
ccm_ag_kkm <- subset(ccm_ag, ccm_ag$pairs=="K-Km")
ccm_ag_kkm$guide <- rep(0, length=nrow(ccm_ag_kkm))
ccm_fig_kkm <- ggplot(ccm_ag_kkm, aes(x=lib_size, y=rho[,3], group=lib_column))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  scale_color_manual(name="Forcing", labels=c("K causes Km", "Km causes K"), values = c("skyblue3", "midnightblue"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  theme(legend.position = c(0.8, 0.12), legend.background = element_rect(fill = "transparent"))
ccm_fig_kkm

#A-K
ccm_ag_kkm <- subset(ccm_ag, ccm_ag$pairs=="K-Km")
ccm_ag_kkm$guide <- rep(0, length=nrow(ccm_ag_kkm))
ccm_fig_kkm <- ggplot(ccm_ag_kkm, aes(x=lib_size, y=rho[,3], group=lib_column))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  scale_color_manual(name="Forcing", labels=c("K causes Km", "Km causes K"), values = c("skyblue3", "midnightblue"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  theme(legend.position = c(0.8, 0.12), legend.background = element_rect(fill = "transparent"))
ccm_fig_kkm


#Km-T
ccm_ag_kmt <- subset(ccm_ag, ccm_ag$pairs=="Km-T")
ccm_ag_kmt$guide <- rep(0, length=nrow(ccm_ag_kmt))
ccm_fig_kmt <- ggplot(ccm_ag_kmt, aes(x=lib_size, y=rho[,3], group=lib_column))+
  #geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  scale_color_manual(name="Forcing", labels=c("T causes Km", "Km causes T"), values = c("grey25","skyblue3"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  theme(legend.position = c(0.8, 0.12), legend.background = element_rect(fill = "transparent"))
ccm_fig_kmt



#**********************
###Time delayed CCM###
#**********************
#Can we make assumptions about indirect interactions/predator-prey delayed effect/syncrony
#(as in Ye et al. 2015 Scientific Reports)?
#Also in the microcosms some species are countable at later times. Could we argue that
#this method helps to deal with that problem?

params   <- expand.grid(lib_column = vars, target_column=vars, tp=-10:10)
params$E <- best_E
params   <- params[params$lib_column != params$target_column,]

ccm_tp_output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i){
  ccm(com_ts, 
      lib = segments, 
      pred = segments, 
      E=params$E[i],  
      lib_column = params$lib_column[i],
      target_column = params$target_column[i], 
      lib_sizes = NROW(data_by_rep[[1]]),
      tp=params$tp[i],
      random_libs = T,
      num_samples = 100, #10000 saved
      silent=T)
}))

#saveRDS(ccm_tp_output, "CCM_delayed_predation.R")
ccm_tp_output <- readRDS("CCM_delayed_predation.R")

ccm_tp_output$direction <- paste(ccm_tp_output$lib_column, "xmap to", ccm_tp_output$target_column)

#CCM_Means
mylist <- split(ccm_tp_output, list(ccm_tp_output$direction, ccm_tp_output$tp))
means <- do.call(rbind, lapply(mylist, function(x) ccm_means(x)))
means$direction <- paste(means$lib_column, "xmap to", means$target_column)

ccm_ag_tp <- with(ccm_tp_output, 
                  aggregate(cbind(rho=rho),
                            list(direction=direction,
                                 lib_size=lib_size,
                                 lib_column=lib_column,
                                 tp=tp),
                            function(x) 
                              quantile(x, c(0.025,pnorm(-1,0,1), 0.5, pnorm(1,0,1), 0.975), 
                                       na.rm=T)))

new_tp_data <- data.frame(
  direction=ccm_ag_tp$direction,
  tp=ccm_ag_tp$tp,
  lower=ccm_ag_tp$rho[,2],
  upper=ccm_ag_tp$rho[,4])

#Here I select the best tp where rho is max. is. 
#In some cases max rho is achieved by chance at higher or lower tps 
#To prevent overfitting I only use the tp where rho max. is when this rho 
#is bigger than the rho at tp=0

#max_tp <- ddply(means, "direction", function(df){
#  cbind(tp=df$tp[which.max(df$rho)],
#        max_rho_min=df$rho[which.max(df$rho)],
#        rho_zero_max=df$rho[which(df$tp==0)],
#        diff=df$rho[which.max(df$rho)]-df$rho[which(df$tp==0)])
#})

#max_tp$best_tp = ifelse(max_tp$diff > 0, max_tp$tp, 0)

#Alternative criteria
#If the lower CI of the max.rho is not overlapping with the rho's higher CI at tp=0 
#use that tp, otherwise 0

max_tp2 <- ddply(new_tp_data, "direction", function(df){
  cbind(tp=df$tp[which.max(df$lower)],
       max_rho_lower=df$lower[which.max(df$lower)],
       rho_zero_upper=df$upper[which(df$tp==0)])
})

max_tp2$best_tp = ifelse(max_tp2$max_rho_lower > max_tp2$rho_zero_upper,
                         max_tp2$tp, 0)


#For plotting
#Dumm
ccm_ag_tp$pairs <-mapvalues(ccm_ag_tp$direction, from = c("A xmap to Am", "Am xmap to A", 
                                                    "A xmap to K", "K xmap to A",
                                                    "A xmap to Km", "Km xmap to A",
                                                    "A xmap to S", "S xmap to A",
                                                    "A xmap to T", "T xmap to A",
                                                    "Am xmap to K", "K xmap to Am",
                                                    "Am xmap to Km", "Km xmap to Am",
                                                    "Am xmap to S", "S xmap to Am", 
                                                    "Am xmap to T", "T xmap to Am",
                                                    "K xmap to Km", "Km xmap to K", 
                                                    "K xmap to S", "S xmap to K",
                                                    "K xmap to T", "T xmap to K", 
                                                    "Km xmap to S", "S xmap to Km", 
                                                    "Km xmap to T", "T xmap to Km",
                                                    "S xmap to T", "T xmap to S"), 
                         to=c("A-Am", "A-Am", "A-K", "A-K", "A-Km", "A-Km", "A-S", "A-S", "A-T", "A-T", 
                              "Am-K", "Am-K", "Am-Km", "Am-Km", "Am-S", "Am-S", "Am-T", "Am-T", 
                              "K-Km","K-Km", "K-S", "K-S", "K-T", "K-T",  "Km-S", "Km-S", "Km-T", 
                              "Km-T", "S-T", "S-T"))
#Plot
as.numeric(as.character(ccm_ag_tp$tp))

time_delay_ccm_fig <- ggplot(ccm_ag_tp, aes(x=tp, y=rho[,3],
                                            color=lib_column))+
  geom_errorbar(aes(x=tp, ymin=rho[,1], ymax=rho[,5]), width=0.5, size=0.4, alpha=0.4)+
  geom_line()+
  facet_wrap(~as.factor(pairs))+
  scale_color_manual(values = c("darkred", "lightpink3", "midnightblue", "skyblue3", 
                                "darkgreen", "grey25"))+
  labs(colour = "Library\ncolumn", x="Cross map lag", y="Cross Map Skill (rho)")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank())

time_delay_ccm_fig

#Example delayed response graph 

#Km-T
ccm_ag_kmt_delayed <- subset(ccm_ag_tp, ccm_ag_tp$pairs=="Km-T")

ccm_tp_kmt <- ggplot(ccm_ag_kmt_delayed, aes(x=tp, y=rho[,3],
                                             color=lib_column))+
  geom_errorbar(aes(x=tp, ymin=rho[,1], ymax=rho[,5]), width=0.5, size=0.4, alpha=0.4)+
  geom_line()+
  scale_color_manual(name="Forcing", labels=c("T causes Km", "Km causes T"), values = c("grey25","skyblue3"))+
  labs(y="Cross Map Skill (rho)", x="Cross map lag")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank())+ 
  theme(legend.position = c(0.8, 0.12), legend.background = element_rect(fill = "transparent"))
ccm_tp_kmt


#**********************
#STOP HERE! #THIS IS PROBABLY SO WRONG!
#CCM E and tp combined
#**********************

params3 <- expand.grid(lib_column = vars, target_column=vars)
params3$E <- best_E
params3 <- params3[params3$lib_column!= params3$target_column,]
params3.ord <-  params3[order(params3$lib_column, params3$target_column),]
params3.ord$tp <- max_tp2$best_tp

ccm_output_best_tp <- do.call(rbind, lapply(seq_len(NROW(params3.ord)), function(i){
  ccm(com_ts, 
      lib = segments, 
      pred = segments, 
      E=params3.ord$E[i],  
      tp=params3.ord$tp[i],
      lib_column = params3.ord$lib_column[i],
      target_column = params3.ord$target_column[i],
      random_libs = T,
      num_samples = 100, #10000 saved
      silent=T)
}))

#saveRDS(ccm_output_best_tp, "CCM_E_tp_combined_predation.R")
ccm_output_best_tp <- readRDS("CCM_E_tp_combined_predation.R")

ccm_output_best_tp$direction <- paste(ccm_output_best_tp$lib_column, "xmap to", 
                                      ccm_output_best_tp$target_column)

ccm_ag_best_tp <- with(ccm_output_best_tp, 
                       aggregate(cbind(rho=rho),
                                 list(direction=direction,
                                      lib_column=lib_column,
                                      lib_size=lib_size),
                                 function(x) quantile(x, c(0.025,
                                                           pnorm(-1,0,1), 0.5,
                                                           pnorm(1,0,1), 0.975),
                                                      na.rm=T)))

ccm_ag_best_tp$pairs <-mapvalues(ccm_ag_best_tp$direction, from = c("A xmap to Am", "Am xmap to A", 
                                                          "A xmap to K", "K xmap to A",
                                                          "A xmap to Km", "Km xmap to A",
                                                          "A xmap to S", "S xmap to A",
                                                          "A xmap to T", "T xmap to A",
                                                          "Am xmap to K", "K xmap to Am",
                                                          "Am xmap to Km", "Km xmap to Am",
                                                          "Am xmap to S", "S xmap to Am", 
                                                          "Am xmap to T", "T xmap to Am",
                                                          "K xmap to Km", "Km xmap to K", 
                                                          "K xmap to S", "S xmap to K",
                                                          "K xmap to T", "T xmap to K", 
                                                          "Km xmap to S", "S xmap to Km", 
                                                          "Km xmap to T", "T xmap to Km",
                                                          "S xmap to T", "T xmap to S"), 
                            to=c("A-Am", "A-Am", "A-K", "A-K", "A-Km", "A-Km", "A-S", "A-S", "A-T", "A-T", 
                                 "Am-K", "Am-K", "Am-Km", "Am-Km", "Am-S", "Am-S", "Am-T", "Am-T", 
                                 "K-Km","K-Km", "K-S", "K-S", "K-T", "K-T",  "Km-S", "Km-S", "Km-T", 
                                 "Km-T", "S-T", "S-T"))


#Same significance criteria as above
signif_tp <- ddply(ccm_output_best_tp,"direction",func1)

signif_tp$V2 <- ifelse(ccm_ag_best_tp[ccm_ag_best_tp$lib_size==max(ccm_ag_best_tp$lib_size),]$rho[,1]>0 &
                         (ccm_ag_best_tp[ccm_ag_best_tp$lib_size==max(ccm_ag_best_tp$lib_size),]$rho[,3]-
                            ccm_ag_best_tp[ccm_ag_best_tp$lib_size==min(ccm_ag_best_tp$lib_size),]$rho[,3])>0,
                       "1","0")

signif_tp

#Plot
mylist2 <- split(ccm_output_best_tp, list(ccm_output_best_tp$direction))
means2 <- do.call(rbind, lapply(mylist2, function(x) ccm_means(x)))
means2$direction <- paste(means2$lib_column, "xmap to", means2$target_column)
means2$guide <- rep(0, length=nrow(means2))


ccm_ag_best_tp$guide <- rep(0, length=nrow(ccm_ag))
ccm_fig_tp <- ggplot(ccm_ag_best_tp, aes(x=lib_size, y=rho[,3], group=lib_column))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  facet_wrap(~pairs)+
  scale_color_manual(name="Library\ncolumn", values = c("darkred", "lightpink3", "midnightblue",   "skyblue3", 
                                "darkgreen", "grey25"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  mytheme+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())
ccm_fig_tp

#********************************************************
###Interaction streght using S-map (Deyle et al. 2016)###
#********************************************************

#Here I created a block data with lags
#If species x caused by y and z and its best E is 4, 
#I create a matrix with x,y,z and x-1 to have the same amount of dimension
#Multivariate predictions
block_data <- make_block(com_ts, t=predation$day, max_lag = 3, lib=segments)
str(block_data)

#Explore the best weighting parameter (nonlinear parameter = theta)
#Best theta is selected based on rho
#Find best theta for all time series 

vars <- c("A","Am","K","Km","S","T")

#what shoud the max theta be?


test_theta_all <- lapply(block_data[,vars], function(vars){
  block_lnlp(vars,
             method = "s-map",
             lib = segments, 
             pred = segments,
             num_neighbors = 0, #use any value < 1 for s-map
             theta = seq(0,4,0.1))
})

best_theta <- sapply(test_theta_all, function(df){
  df$theta[which.max(round(df$rho,2))]
})
best_theta


par(mfrow=c(3,2))
for(target in names(test_theta_all)){
  plot(test_theta_all[[target]]$theta, test_theta_all[[target]]$rho, type="l", 
       xlab="theta", ylab="rho", main=target)
}

#Do S-map analysis with the best theta for A
#A E=5

#I used all the variables but I guess for the manuscript
#we have to construct the interaction matrix according to the 
#ccm results. 
#If max. E is 5 and variable forced by only 3 other variables, 
#remaining 2 dimentions should be filled with the lags of the variables

#Maybe this can be written as a function. CCM 1/0 matrix, E and block data... 

smap_res_A <-  block_lnlp(block_data[,vars],
                        method = "s-map",
                        num_neighbors = 0, 
                        lib= segments, 
                        pred=segments,
                        theta = best_theta["A"],
                        #columns = vars,
                        target_column = "A",
                        silent = T,
                        save_smap_coefficients = T, 
                        first_column_time = F) 


#Time series of fluctuating interaction strength
smap_coef_A <- as.data.frame(smap_res_A$smap_coefficients[[1]])

#colnames(smap_coef) ???
#These colnames are confusing!!!
#I took the mean of the predictions for 3 replicates 
smap_coef_A_mean <- data.frame((smap_coef_A[1:31,]+smap_coef_A[33:63,]+smap_coef_A[65:95,])/3)

#Smap for Am
#Am E=4
smap_res_Am <-  block_lnlp(block_data[,vars],
                           method = "s-map",
                           num_neighbors = 0, 
                           lib= segments, 
                           pred=segments,
                           theta = best_theta["Am"],
                           #columns = c("Am", "K", "Km", "S"),
                           target_column = "Am",
                           silent = T,
                           save_smap_coefficients = T, 
                           first_column_time = F) 

#Time series of fluctuating interaction strength
smap_coef_Am <- as.data.frame(smap_res_Am$smap_coefficients[[1]])

smap_coef_Am_mean <- data.frame((smap_coef_Am[1:31,]+smap_coef_Am[33:63,]+smap_coef_Am[65:95,])/3)

#Smap for K
#None causes K

smap_res_K <-  block_lnlp(block_data[,vars],
                           method = "s-map",
                           num_neighbors = 0, 
                           lib= segments, 
                           pred=segments,
                           theta = best_theta["K"],
                           #columns = 
                           target_column = "K",
                           silent = T,
                           save_smap_coefficients = T, 
                           first_column_time = F) 

#Time series of fluctuating interaction strength
smap_coef_K <- as.data.frame(smap_res_K$smap_coefficients[[1]])

smap_coef_K_mean <- data.frame((smap_coef_K[1:31,]+smap_coef_K[33:63,]+smap_coef_K[65:95,])/3)


#Smap for Km
#Km E=4

smap_res_Km <-  block_lnlp(block_data[,vars],
                           method = "s-map",
                           num_neighbors = 0, 
                           lib= segments, 
                           pred=segments,
                           theta = best_theta["Km"],
                           #columns = vars,
                           target_column = "Km",
                           silent = T,
                           save_smap_coefficients = T, 
                           first_column_time = F) 

#Time series of fluctuating interaction strength
smap_coef_Km <- as.data.frame(smap_res_Km$smap_coefficients[[1]])

smap_coef_Km_mean <- data.frame((smap_coef_Km[1:31,]+smap_coef_Km[33:63,]+smap_coef_Km[65:95,])/3)

#Smap for S
#S E=4

smap_res_S <-  block_lnlp(block_data[,vars],
                          method = "s-map",
                          num_neighbors = 0, 
                          lib= segments, 
                          pred=segments,
                          theta = best_theta["S"],
                          #columns= vars,
                          target_column = "S",
                          silent = T,
                          save_smap_coefficients = T, 
                          first_column_time = F) # save S-map coefficients

#Time series of fluctuating interaction strength
smap_coef_S <- as.data.frame(smap_res_S$smap_coefficients[[1]])

smap_coef_S_mean<-data.frame((smap_coef_S[1:31,]+smap_coef_S[33:63,]+smap_coef_S[65:95,])/3)

#T E=2

#Smap for T
smap_res_T <-  block_lnlp(block_data[,vars],
                        method = "s-map",
                          num_neighbors = 0,
                          lib= segments, 
                          pred=segments,
                          theta = best_theta["T"],
                          #columns =vars, 
                          target_column = "T",
                          silent = T,
                          save_smap_coefficients = T, 
                          first_column_time = F) 

#Time series of fluctuating interaction strength
smap_coef_T <- as.data.frame(smap_res_T$smap_coefficients[[1]])

smap_coef_T_mean <-data.frame((smap_coef_T[1:31,]+smap_coef_T[33:63,]+smap_coef_T[65:95,])/3)

# Here I calculate the mean interaction strenght. I am not sure if I have to standardize between -1-1. 
#And I am also not sure about the columns are in the original order of the species. 
#So I assume that each column is the partial derivatives of the target species divided by another species (effect of a species on the target species). 

names(smap_coef_A_mean)   <- c("A", "Am", "K", "Km", "S", "T")
names(smap_coef_Am_mean)  <- c("A", "Am", "K", "Km", "S", "T")
names(smap_coef_K_mean)   <- c("A", "Am", "K", "Km", "S", "T")
names(smap_coef_Km_mean)  <- c("A", "Am", "K", "Km", "S", "T")
names(smap_coef_S_mean)   <- c("A", "Am", "K", "Km", "S", "T")
names(smap_coef_T_mean)   <- c("A", "Am", "K", "Km", "S", "T")

all_coef <- rbind(
  melt(cbind(smap_coef_A_mean[,1:6], target=rep("A", nrow(smap_coef_A_mean)))),
  melt(cbind(smap_coef_Am_mean[,1:6], target=rep("Am", nrow(smap_coef_Am_mean)))),
  melt(cbind(smap_coef_K_mean[,1:6], target=rep("K", nrow(smap_coef_K_mean)))),
  melt(cbind(smap_coef_Km_mean[,1:6], target=rep("Km", nrow(smap_coef_Km_mean)))),
  melt(cbind(smap_coef_S_mean[,1:6], target=rep("S", nrow(smap_coef_S_mean)))),
  melt(cbind(smap_coef_T_mean[,1:6], target=rep("T", nrow(smap_coef_T_mean))))
)


#Plot 
palette2   <- c("darkred", "lightpink3", "midnightblue", "skyblue3", "darkgreen", "grey25")

ggplot(all_coef,aes(y=value, x=as.factor(as.character(target)), color=target))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_jitter(alpha=0.5)+
  geom_boxplot()+
  facet_wrap(.~variable, scales="free", drop = T)+
  scale_color_manual(values=palette2)+
  labs(x="Effect of", y="Interaction strength")+
  mytheme+
  theme(legend.position="none")+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Plot 
all_coef$time <- rep(3:33, times=36)

ggplot(all_coef,aes(y=value, x=as.numeric(time), group=target, color=target))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line()+
  facet_wrap(.~variable, scales="free_y", drop = T, strip.position = "top")+
  scale_color_manual(name="Effect of", values=palette2)+
  labs(x="Time", y="Interaction strength")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


aggdata <- aggregate(all_coef$value, by = list(target=all_coef$target, 
                                               variable=all_coef$variable), FUN = mean, na.rm=T)


#Here is an overview of the positive-negative interactions
aggdata$sign <- ifelse(aggdata$x>0, "positive",
                       ifelse(aggdata$x<0, "negative", "none"))


sign <- ggplot(aggdata, aes(y = variable, x = target, fill = sign)) +
  geom_tile(color="white")+
  labs(color="Interaction", y="Target", x="Effect", title="Type of interactions")+
  scale_fill_manual(name="Interaction\ntype", values=c("#999999", "#E69F00"))+
  mytheme+
  theme(legend.position = "none")

 freq <- ggplot(aggdata, aes(x=abs(x),fill=sign)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  scale_color_manual(name= "Type", values=c("#999999", "#E69F00"))+
  scale_fill_manual(name= "Type", values=c("#999999", "#E69F00"))+
  labs(title="Skewness",x="Interaction Strenght", y = "Frequency")+
  mytheme+
  theme(plot.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
   theme(legend.position = c(0.8, 0.8))


 plot_grid(sign, freq)
 

#Example interaction plot 
coef_example <- subset(all_coef, all_coef$variable =='A' | all_coef$variable =='K' | all_coef$variable =='T')

ggplot(coef_example,aes(y=value, x=as.numeric(time), group=target, color=target))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line()+
  facet_wrap(.~variable, scales="free_y", drop = T, strip.position = "top")+
  scale_color_manual(name="Effect of", values=palette2)+
  labs(x="Time", y="Interaction strength")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#Example interaction plot 
coef_example <- subset(all_coef, all_coef$variable =='Km' | all_coef$variable =='T' )

ggplot(coef_example,aes(y=value, x=as.numeric(time), group=target, color=target))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(size=1)+
  facet_wrap(.~variable, scales="free_y", drop = T, strip.position = "top")+
  scale_color_manual(name="Effect of", values=palette2)+
  labs(x="Time", y="Interaction strength")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#********************
# Dynamic stability 
#********************
#Interaction matrix for each time point 
#This is so wrong 
#I used all the interactions for now 

#But interaction strenght of the ones CCM 0 should be fixed to 0 
#removing doesn't work, because the matrix should be quadratic 
#if there are only a couple of variables it doesn't work 

stab_data <- as.data.frame(rbind(
  cbind(time = seq(3,33), smap_coef_A_mean[,c(1,2,3,4,5,6)], target=rep("A", nrow(smap_coef_A_mean))),
  cbind(time = seq(3,33), smap_coef_Am_mean[,c(1,2,3,4,5,6)], target=rep("Am", nrow(smap_coef_Am_mean))),
  cbind(time = seq(3,33), smap_coef_K_mean[,c(1,2,3,4,5,6)], target=rep("K", nrow(smap_coef_K_mean))),
  cbind(time = seq(3,33), smap_coef_Km_mean[,c(1,2,3,4,5,6)], target=rep("Km", nrow(smap_coef_Km_mean))),
  cbind(time = seq(3,33), smap_coef_S_mean[,c(1,2,3,4,5,6)], target=rep("S", nrow(smap_coef_S_mean))),
  cbind(time = seq(3,33), smap_coef_T_mean[,c(1,2,3,4,5,6)], target=rep("T", nrow(smap_coef_T_mean)))
))

#Order and split data into time points
stab_data_sort <- stab_data[order(stab_data$time),]
stab_split <- split(stab_data_sort[,2:7], stab_data_sort$time)

#Calculating eigenvalues
eigenvalues <- sapply(stab_split, function(x){
  max(Re(eigen(x)$values))
}) #it doesn't matter if you transfom the matrix or not 

stability <- melt(eigenvalues)
stability$time <- rownames(stability)

#Interaction indexes
time_mean <- sapply(stab_split, function(x){mean(abs(x[row(x)!=col(x)]))})
time_median <- sapply(stab_split, function(x){median(abs(x[row(x)!=col(x)]))})
time_max <- sapply(stab_split, function(x){max(abs(x[row(x)!=col(x)]))})

#Weak interactions index
weak <- as.vector(time_median)/as.vector(time_max) 

mean.interactions <- melt(time_mean)
mean.interactions$time <- rownames(mean.interactions)

weak.interactions <- melt(weak)
weak.interactions$time <- rownames(stability)

#Diversity index (exponential Shannon)
pred_split   <- split(predation, predation$replicate)
pred_shannon <- apply(predation[,3:8],1, function(x) {exp(diversity(x, index ="shannon"))})
pred_shannon_data <- cbind(predation[,1:2], pred_shannon)

library(dplyr)
grouped_shannon <- group_by(pred_shannon_data, day) %>% summarise(mean=mean(pred_shannon))

#Coefficient of variation 
cv <- function(x) {sd(x)/mean(x)}
grouped_abun  <- group_by(predation_ts, day, variable) %>% summarise(cv=cv(value+0.1)) %>%
  summarise(mean_cv=mean(cv))

#Plots 
stab <- ggplot(stability,aes(y=value, x=as.numeric(time)))+
  geom_hline(yintercept = 1, linetype="dashed", color="red",size=1.5)+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y="Stability")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

mean.int <- ggplot(mean.interactions,aes(y=value, x=as.numeric(time)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y="Mean interaction\nstrength")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

weak.int <- ggplot(weak.interactions,aes(y=value, x=as.numeric(time)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  geom_segment(mapping=aes(x=3, y=0.25, xend=3, yend=0.18), 
               size=0.7, color="black", 
               arrow = arrow(length = unit(0.4, "cm")))+
  annotate("text", label = "Dominance of\nweak interactions", x = 4.5, y = 0.22, 
           color = "black", angle=90)+
  labs(x="Time", y="Median:maximum\ninteraction strenght")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pred_cv <- ggplot(grouped_abun,aes(y=mean_cv, x=as.numeric(day)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y="Coefficient of\nvariation")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pred_diversity <- ggplot(grouped_shannon,aes(y=mean, x=as.numeric(day)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y="Diversity index\nHill 1")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Abundances
abun.pred  <- as.data.frame(subset(dfc_pred, dfc_pred$variable=="T"))
abun.morph <- as.data.frame(subset(dfc_pred, dfc_pred$variable=="Km"))

#This is a bit too much! 
#Anyways here's wonderwall 

#See what causes dynamic stability
all_indexes <- as.data.frame(cbind(dyn.stab=stability$value, mean.int=mean.interactions$value, 
                     weak.int=weak.interactions$value, div=grouped_shannon$mean, var=grouped_abun$mean_cv, 
                     pred.abun=abun.pred$value, morph.abun=abun.morph$value))

all_indexes_norm <- as.data.frame(apply(all_indexes, 2, function(x) (x - min(x))/(max(x)-min(x))))

#Find best E
simplex_stab <- simplex(all_indexes_norm$dyn.stab, E=1:6, tau=1)

par(mfrow=c(1,1))
plot(simplex_stab$E, simplex_stab$rho, type="l", 
       xlab="Embedding Dimension (E)", ylab="Forecast Skill (rho)")
#E=1

#*********
###CCM###
#********

#CCM for tp=0 which is default
#Effect of mean interactions

ccm_output_stab_int <- ccm(all_indexes_norm, 
                            E=1,  
                            lib_column = "dyn.stab",
                            target_column = "mean.int",
                            random_libs = T,
                            num_samples = 10000,
                            silent=T)


ccm_output_stab_int$direction <- paste(ccm_output_stab_int$lib_column, "xmap to", ccm_output_stab_int$target_column)

ccm_ag_stab_int <- with(ccm_output_stab_int, aggregate(cbind(rho=rho),
                                     list(direction=direction,
                                          lib_column=lib_column,
                                          lib_size=lib_size),
                                     function(x) quantile(x, c(0.025,
                                                               pnorm(-1,0,1), 0.5, 
                                                               pnorm(1,0,1), 0.975), 
                                                          na.rm=T)))

#Plot 
ccm_ag_stab_int$guide <- rep(0, length=nrow(ccm_ag_stab_int))
ccm_fig_stab_int <- ggplot(ccm_ag_stab_int, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Mean interaction strength")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_int


#Effect of weak interactions
ccm_output_stab_weak <- ccm(all_indexes_norm, 
                           E=1,  
                           lib_column = "dyn.stab",
                           target_column = "weak.int",
                           random_libs = T,
                           num_samples = 10000,
                           silent=T)


ccm_output_stab_weak$direction <- paste(ccm_output_stab_weak$lib_column, "xmap to", ccm_output_stab_weak$target_column)

ccm_ag_stab_weak <- with(ccm_output_stab_weak, aggregate(cbind(rho=rho),
                                                       list(direction=direction,
                                                            lib_column=lib_column,
                                                            lib_size=lib_size),
                                                       function(x) quantile(x, c(0.025,
                                                                                 pnorm(-1,0,1), 0.5, 
                                                                                 pnorm(1,0,1), 0.975), 
                                                                            na.rm=T)))

#Plot 
ccm_ag_stab_weak$guide <- rep(0, length=nrow(ccm_ag_stab_weak))
ccm_fig_stab_weak <- ggplot(ccm_ag_stab_weak, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Weak interactions")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_weak


#Effect of diversity
ccm_output_stab_div <- ccm(all_indexes_norm, 
                            E=1,  
                            lib_column = "dyn.stab",
                            target_column = "div",
                            random_libs = T,
                            num_samples = 10000,
                            silent=T)


ccm_output_stab_div$direction <- paste(ccm_output_stab_div$lib_column, "xmap to", ccm_output_stab_div$target_column)

ccm_ag_stab_div <- with(ccm_output_stab_div, aggregate(cbind(rho=rho),
                                                         list(direction=direction,
                                                              lib_column=lib_column,
                                                              lib_size=lib_size),
                                                         function(x) quantile(x, c(0.025,
                                                                                   pnorm(-1,0,1), 0.5, 
                                                                                   pnorm(1,0,1), 0.975), 
                                                                              na.rm=T)))

#Plot 
ccm_ag_stab_div$guide <- rep(0, length=nrow(ccm_ag_stab_div))
ccm_fig_stab_div <- ggplot(ccm_ag_stab_div, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Diversity index")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_div



#Effect of coefficient of variation #no causation
#Effect of morphotype abundance #no causation


#Effect of predator abundance
ccm_output_stab_pred <- ccm(all_indexes_norm, 
                           E=1,  
                           lib_column = "dyn.stab",
                           target_column = "pred.abun",
                           random_libs = T,
                           num_samples = 10000,
                           silent=T)


ccm_output_stab_pred$direction <- paste(ccm_output_stab_pred$lib_column, "xmap to", ccm_output_stab_pred$target_column)

ccm_ag_stab_pred <- with(ccm_output_stab_pred, aggregate(cbind(rho=rho),
                                                       list(direction=direction,
                                                            lib_column=lib_column,
                                                            lib_size=lib_size),
                                                       function(x) quantile(x, c(0.025,
                                                                                 pnorm(-1,0,1), 0.5, 
                                                                                 pnorm(1,0,1), 0.975), 
                                                                            na.rm=T)))

#Plot 
ccm_ag_stab_pred$guide <- rep(0, length=nrow(ccm_ag_stab_pred))
ccm_fig_stab_pred <- ggplot(ccm_ag_stab_pred, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Predator abundance")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_pred


plot_grid(ccm_fig_stab_int, ccm_fig_stab_weak, ccm_fig_stab_div,
          ccm_fig_stab_pred, ncol=2)

#Interaction strenghts 

all_indexes_norm_block_data <- make_block(all_indexes_norm,max_lag = 3)

all_indexes_norm_selected <- all_indexes_norm_block_data[,
                                        c("dyn.stab","mean.int", "weak.int", "div", "pred.abun")]


test_theta_all <-
  block_lnlp(all_indexes_norm_selected$dyn.stab,
             method = "s-map",
             num_neighbors = 0, 
             theta = seq(0,4,0.1))

test_theta_all$theta[which.max(round(test_theta_all$rho,2))]


smap_stab <-  block_lnlp(all_indexes_norm_selected,
                          method = "s-map",
                          num_neighbors = 0,
                          theta = 0,
                          target_column = "dyn.stab",
                          silent = T,
                          save_smap_coefficients = T, 
                          first_column_time = F) 


#Time series of fluctuating interaction strength
smap_stab_dyn <- as.data.frame(smap_stab$smap_coefficients[[1]])

smap_stab_dyn_sub <- melt(smap_stab_dyn[,2:5])

contribution <- ggplot(smap_stab_dyn_sub, aes(x=variable, y=abs(value), fill=variable))+
  geom_boxplot()+
  ylab("Influences on stability")+
  scale_x_discrete(labels = c("Mean interspecific\ninteraction strenght", 'Dominance of\nweak interactions', 'Diversity', 
                              "Predator abundance"))+
  scale_fill_manual(values = c("#999999", "#999999", "#E69F00", "#999999"))+
  mytheme+
  theme(legend.position= "none",
        axis.title.x = element_blank(), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


#Combined Interaction plot 
plot_grid(stab, mean.int, weak.int, pred_cv, pred_diversity, contribution,
          ncol=2)


#Little test 

#Predator-prey ratio 
#This is nice but it is a bit dumm. 
#First I say Km doesn't cause K, and than calculate 
#interaction strengths and show this magic graph
#So wrong...

predation <- as.data.frame(cbind(day=predation.data$day, 
                                 ddply(predation.data[,2:7], 
                                       .(replicate=as.factor(predation.data$replicate)), 
                                       colwise(na.approx))))

prey <- log(predation$K+predation$Km+1)
ratio <-log(predation$T)/prey

#Partial derivatives K/Km
derivs <- smap_coef_K$c_4
plot(derivs~ratio)

test_data <- as.data.frame(cbind(derivs=derivs, pred.abun=predation$T, ratio))

test_fig <- ggplot(test_data, aes(x=pred.abun, y=derivs))+
  geom_hline(yintercept=0, linetype="dashed", color="red", size=1.2)+
  geom_point(shape=8)+
  stat_smooth(method="lm", formula=y~poly(x,2))+
  ylab(expression(atop("Effect of Km on K",paste(competition %<->% mutualism))))+
  xlab("Predator abundance")+
  mytheme+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.background = element_blank()
        )


#Draw the networks with the CCM results/significance 

#Create incidence matrices
#This network drawing thing should be optimized
#the easiest and most customable network ilustration 
#was geom_net of ggplot, but I had to modify in inkScape

#signif
#signif_tp
incidences <- (signif %>% separate(direction, 
                                   c("variable", "rest","rest", "target")))[,c(1,4,6)]

#Here I messed up 
#add the positive-negative information 
#network.data <- (merge(incidences, aggdata, by=c("variable", "target")))[,c(1,2,3,5)]
#names(network.data) <- c("to", "from", "significance", "type")

#Subset the significant causations for plotting purposes 
#network.data.subset <- subset(network.data, network.data$significance==1)

#add the positive-negative information 
network.data <- merge(incidences, aggdata, by=c("variable", "target"), all=T)
#names(network.data) <- c("to", "from", "significance", "strength","type")

#Subset the significant causations for plotting purposes 
#network.data.subset2_c <- subset(network.data2_c, network.data2_c$significance==1)

network.data$variable <-mapvalues(network.data$variable, from = c("A", "Am","K", "Km","S","T"), 
                                 to=c("6A", "2Am", "5K", "3Km", "4S", "1T"))
network.data$target <-mapvalues(network.data$target, from = c("A", "Am","K", "Km","S","T"), 
                                  to=c("6A", "2Am", "5K", "3Km", "4S", "1T"))
palette5   <- c("grey25", "lightpink3", "skyblue3", "darkgreen", "midnightblue", "darkred")

net <- ggplot(data = network.data, aes(from_id = target, to_id = variable)) +
  geom_net(layout.alg = 'circle', labelon = TRUE, 
           size = 11, directed = TRUE, vjust = 0.5, ealpha = 0.7,
           arrowsize = 1, arrowgap = 0.07, singletons= TRUE, labelcolour = "white", 
           fontsize = 5, curvature=0.1,
           aes(color=target, shape=target,linewidth = rescale(abs(x), c(0.4,6))))+ 
  theme_net()+
  scale_colour_manual(values=palette5)+
  scale_linetype_manual(values=c("twodash", "solid"))+
  scale_shape_manual(values=c(17,16,16,16,16,16))+
  theme(panel.background = element_rect(color = "grey"))+
  theme(legend.position = "none")

#Rest is the copy of the code above. 

#*********************************************
###SAME FOR COMMUNUNITIES WITHOUT PREDATOR###
#*********************************************

#Separate time column from variables 
#Split data by replicates
#I use all the segments for prediction 

vars_c <- c("A","Am","S")

com_ts_c <- competition[,vars_c]
data_by_rep_c <- split(com_ts_c, competition$replicate)
segments_end_c <- cumsum(sapply(data_by_rep_c, NROW))
segments_begin_c <- c(1, segments_end_c[-length(segments_end_c)]+1)
segments_c <- cbind(segments_begin_c,segments_end_c) #Use all segments for predictions

#Predictability & nonlinearity
#Find best E
simplex_out_c <- lapply(vars_c, function(vars_c){
  simplex(competition[,c("day", vars_c)], E=1:3, lib=segments_c, pred=segments_c)
})

names(simplex_out_c) <- vars_c

par(mfrow=c(3,2))
for(var in names(simplex_out_c)){
  plot(simplex_out_c[[var]]$E, simplex_out_c[[var]]$rho, type="l", 
       xlab="Embedding Dimension (E)", ylab="Forecast Skill (rho)", main = var)
}

max_E_c <- sapply(simplex_out_c, function(df){
  df$E[which.max(round(df$rho,1))]
})
max_E_c

#Choose the smaller E if the difference is minor 
best_E_c <- as.integer(c(1,2,3))
names(best_E_c) <- vars_c


##Nonlinearity
smap_out_c <- lapply(vars_c, function(var){
  s_map(competition[,c("day", var)], E=best_E[var], lib= segments)
})
names(smap_out_c) <- names(simplex_out_c)

#par(mfrow=c(2,2))
for(var in names(smap_out_c)){
  plot(smap_out_c[[var]]$theta, smap_out_c[[var]]$rho, type="l", 
       xlab="Nonlinearity (theta)", ylab="Forcast Skill (rho)", main=var)
}

#*********
###CCM###
#*********

#CCM for tp=0 which is default
vars_c <- c("A","Am","S")
params2_c <- expand.grid(lib_column = vars_c, target_column=vars_c)
params2_c$E <- best_E_c
params2_c <- params2_c[params2_c$lib_column!= params2_c$target_column,]
params2.ord_c <-  params2_c[order(params2_c$lib_column, params2_c$target_column),]

ccm_output_c <- do.call(rbind, lapply(seq_len(NROW(params2.ord_c)), function(i){
  ccm(com_ts_c, 
      lib = segments_c, 
      pred = segments_c, 
      E=params2.ord_c$E[i],  
      lib_column = params2.ord_c$lib_column[i],
      target_column = params2.ord_c$target_column[i],
      random_libs = T,
      num_samples = 100, # Saved 10000
      silent=T)
}))

#saveRDS(ccm_output_c, "CCM_tp0_competition.R")
ccm_output_c <- readRDS("CCM_tp0_competition.R")


ccm_output_c$direction <- paste(ccm_output_c$lib_column, "xmap to", ccm_output_c$target_column)

ccm_ag_c <- with(ccm_output_c, aggregate(cbind(rho=rho),
                                     list(direction=direction,
                                          lib_size=lib_size, 
                                          lib_column=lib_column),
                                     function(x) quantile(x, c(0.025,
                                                               pnorm(-1,0,1), 0.5, 
                                                               pnorm(1,0,1), 0.975), na.rm=T)))


ccm_ag_c$pairs <-mapvalues(ccm_ag_c$direction, from = c("A xmap to Am", "Am xmap to A", 
                                                                    "A xmap to S", "S xmap to A",
                                                                    "Am xmap to S", "S xmap to Am"), 
                                 to=c("A-Am", "A-Am", "A-S", "A-S", "Am-S", "Am-S"))

#Significance criteria according to Clark et al. Ecology paper 
func1 <- function(data){1-mean((data[data$lib_size==min(data$lib_size),]$rho<data[data$lib_size==max(data$lib_size),]$rho)
                               &(data[data$lib_size==max(data$lib_size),]$rho>0), na.rm=T)
}
signif_c <- ddply(ccm_output_c,"direction",func1)

#Significance criteria, less conservative
#If terminal rho is bigger than 0 and the difference between terminal rho and initial rho is
#bigger than 0 (convergence)

signif_c$V2 <- ifelse(ccm_ag_c[ccm_ag_c$lib_size==max(ccm_ag_c$lib_size),]$rho[,1]>0 &
                      (ccm_ag_c[ccm_ag_c$lib_size==max(ccm_ag_c$lib_size),]$rho[,1]-
                         ccm_ag_c[ccm_ag_c$lib_size==min(ccm_ag_c$lib_size),]$rho[,1])>0, 
                      "1","0")
signif_c

#Plot 
ccm_ag_c$guide <- rep(0, length=nrow(ccm_ag_c))

ccm_fig_c <- ggplot(ccm_ag_c, aes(x=lib_size, y=rho[,3], group=lib_column))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  facet_wrap(~pairs)+
  scale_color_manual(name= "Library\ncolumn", values = c("darkred", "lightpink3", "darkgreen"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_c

#*********************
###Time delayed CCM###
#*********************

params_c   <- expand.grid(lib_column = vars_c, target_column=vars_c, tp=-10:10)
params_c$E <- best_E_c
params_c   <- params_c[params_c$lib_column!= params_c$target_column,]

ccm_tp_output_c <- do.call(rbind, lapply(seq_len(NROW(params_c)), function(i){
  ccm(com_ts_c, 
      lib = segments_c, 
      pred = segments_c, 
      E=params_c$E[i],  
      lib_column = params_c$lib_column[i],
      target_column = params_c$target_column[i], 
      lib_sizes = NROW(data_by_rep_c[[1]]),
      tp=params_c$tp[i],
      random_libs = T,
      num_samples = 100, #10000
      silent=T)
}))


#saveRDS(ccm_tp_output_c, "CCM_delayed_competition.R")
ccm_tp_output_c <- readRDS("CCM_delayed_competition.R")

ccm_tp_output_c$direction <- paste(ccm_tp_output_c$lib_column, "xmap to", ccm_tp_output_c$target_column)

#CCM_Means
mylist_c <- split(ccm_tp_output_c, list(ccm_tp_output_c$direction, ccm_tp_output_c$tp))
means_c <- do.call(rbind, lapply(mylist_c, function(x) ccm_means(x)))
means_c$direction <- paste(means_c$lib_column, "xmap to", means_c$target_column)

# Here I select the best tp where rho is max. is. 
# In some cases max rho is achieved by chance at higher or lower tps 
# To prevent overfitting I only use the tp where rho max. is when this rho 
# is bigger than the rho at tp=0

#max_tp_c <- ddply(means_c, "direction", function(df){
#  cbind(tp=df$tp[which.max(df$rho)],
#        max_rho_min=df$rho[which.max(df$rho)],
#        rho_zero_max=df$rho[which(df$tp==0)],
#        diff=df$rho[which.max(df$rho)]-df$rho[which(df$tp==0)])
#})

#max_tp_c$best_tp_c = ifelse(max_tp_c$diff > 0, max_tp_c$tp, 0)

ccm_ag_tp_c <- with(ccm_tp_output_c, 
                  aggregate(cbind(rho=rho),
                            list(direction=direction,
                                 lib_size=lib_size,
                                 lib_column=lib_column,
                                 tp=tp),
                            function(x) 
                              quantile(x, c(0.025,pnorm(-1,0,1), 0.5, pnorm(1,0,1), 0.975), 
                                       na.rm=T)))

new_tp_data_c <- data.frame(
  direction=ccm_ag_tp_c$direction,
  tp=ccm_ag_tp_c$tp,
  lower=ccm_ag_tp_c$rho[,2],
  upper=ccm_ag_tp_c$rho[,4])

#Alternative criteria
#If the lower CI of the max.rho is not overlapping with the rho's higher CI at tp=0 
#use that tp, otherwise 0

max_tp2_c <- ddply(new_tp_data_c, "direction", function(df){
  cbind(tp=df$tp[which.max(df$lower)],
        max_rho_lower=df$lower[which.max(df$lower)],
        rho_zero_upper=df$upper[which(df$tp==0)])
})

max_tp2_c$best_tp = ifelse(max_tp2_c$max_rho_lower > max_tp2_c$rho_zero_upper, 
                           max_tp2_c$tp, 0)

#For plotting
ccm_ag_tp_c$pairs <- mapvalues(ccm_ag_tp_c$direction, from = c("A xmap to Am", "Am xmap to A", 
                                                            "A xmap to S", "S xmap to A",
                                                            "Am xmap to S", "S xmap to Am"), 
                               to=c("A-Am", "A-Am", "A-S", "A-S", "Am-S", "Am-S"))


#Plot
as.numeric(as.character(ccm_ag_tp_c$tp))

time_delay_ccm_fig_c <- ggplot(ccm_ag_tp_c, aes(x=tp, y=rho[,3],
                                            color=lib_column))+
  geom_errorbar(aes(x=tp, ymin=rho[,1], ymax=rho[,5]), width=0.5, size=0.4, alpha=0.4)+
  geom_line()+
  facet_wrap(~as.factor(pairs))+
  scale_color_manual(values = c("darkred", "lightpink3", "darkgreen"))+
  labs(colour = "Library\ncolumn", x="Cross map lag", y="Cross Map Skill (rho)")+
  mytheme+ 
  theme(plot.background = element_blank(),
                panel.grid.major = element_blank())

time_delay_ccm_fig_c

#**********************
#STOP HERE!
#CCM E and tp combined
#**********************

params3_c <- expand.grid(lib_column = vars_c, target_column=vars_c)
params3_c$E <- best_E_c
params3_c <- params3_c[params3_c$lib_column!= params3_c$target_column,]
params3.ord_c <-  params3_c[order(params3_c$lib_column, params3_c$target_column),]
params3.ord_c$tp <- max_tp2_c$best_tp

ccm_output_best_tp_c <- do.call(rbind, lapply(seq_len(NROW(params3.ord_c)), function(i){
  ccm(com_ts_c, 
      lib = segments_c, 
      pred = segments_c, 
      E=params3.ord_c$E[i],  
      tp=params3.ord_c$tp[i],
      lib_column = params3.ord_c$lib_column[i],
      target_column = params3.ord_c$target_column[i],
      random_libs = T,
      num_samples = 10000, #Saved 10000
      silent=T)
}))

saveRDS(ccm_output_best_tp_c, "CCM_E_tp_combined_competition.R")
#ccm_output_best_tp_c <- readRDS("CCM_E_tp_combined_competition.R")

ccm_output_best_tp_c$direction <- paste(ccm_output_best_tp_c$lib_column, "xmap to", 
                                      ccm_output_best_tp_c$target_column)


ccm_ag_best_tp_c <- with(ccm_output_best_tp_c, 
                       aggregate(cbind(rho=rho),
                                 list(direction=direction,
                                      lib_size=lib_size, 
                                      lib_column=lib_column),
                                 function(x) quantile(x, c(0.025,
                                                           pnorm(-1,0,1), 0.5,
                                                           pnorm(1,0,1), 0.975),
                                                      na.rm=T)))

ccm_ag_best_tp_c$pairs <- mapvalues(ccm_ag_best_tp_c$direction, from = c("A xmap to Am", "Am xmap to A", 
                                                                   "A xmap to S", "S xmap to A",
                                                                   "Am xmap to S", "S xmap to Am"), 
                                 to=c("A-Am", "A-Am", "A-S", "A-S", "Am-S", "Am-S"))

#Same significance criteria as above
signif_tp_c <- ddply(ccm_output_best_tp_c,"direction",func1)

signif_tp_c$V2 <- ifelse(ccm_ag_best_tp_c[ccm_ag_best_tp_c$lib_size==max(ccm_ag_best_tp_c$lib_size),]$rho[,1]>0 &
                         (ccm_ag_best_tp_c[ccm_ag_best_tp_c$lib_size==max(ccm_ag_best_tp_c$lib_size),]$rho[,3]-
                            ccm_ag_best_tp_c[ccm_ag_best_tp_c$lib_size==min(ccm_ag_best_tp_c$lib_size),]$rho[,3])>0,
                       "1","0")

signif_tp_c

#Plot

ccm_ag_best_tp_c$guide <- rep(0, length=nrow(ccm_ag_best_tp_c))

ccm_fig_tp_c <- ggplot(ccm_ag_best_tp_c, aes(x=lib_size, y=rho[,3], group=lib_column))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(aes(color=lib_column))+
  facet_wrap(~pairs)+
  scale_color_manual(name= "Library\ncolumn", values = c("darkred", "lightpink3", "darkgreen"))+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_tp_c

#********************************************************
###Interaction streght using S-map (Deyle et al. 2016)###
#********************************************************

#Here I created a block data using max. one lag 

#Multivariate predictions
block_data_c <- make_block(com_ts_c, t=competition$day, max_lag = 3, lib=segments_c)
str(block_data_c)

#Explore the best weighting parameter (nonlinear parameter = theta)
#Best theta is selected based on rho
#Find best theta for all time series 

vars_c <- c("A","Am","S")

test_theta_all_c <- lapply(block_data_c[,vars_c], function(vars_c){
  block_lnlp(vars_c,
             method = "s-map",
             lib = segments_c, 
             pred = segments_c,
             num_neighbors = 0, 
             theta = seq(0,4,0.1))
})

best_theta_c <- sapply(test_theta_all_c, function(df){
  df$theta[which.max(round(df$rho,2))]
})
best_theta_c

par(mfrow=c(1,3))
for(target in names(test_theta_all_c)){
  plot(test_theta_all_c[[target]]$theta, test_theta_all_c[[target]]$rho, type="l", 
       xlab="theta", ylab="rho", main=target)
}

#Do S-map analysis with the best theta for A
#A E=1

smap_res_A_c <-  block_lnlp(block_data_c[,vars_c],
                        method = "s-map",
                        num_neighbors = 0,
                        lib= segments_c, 
                        pred=segments_c,
                        theta = best_theta_c[1],
                        #columns = vars_c,
                        target_column = "A",
                        silent = T,
                        save_smap_coefficients = T, 
                        first_column_time = F) 


#Time series of fluctuating interaction strength
smap_coef_A_c <- as.data.frame(smap_res_A_c$smap_coefficients[[1]])
#colnames(smap_coef) ???

smap_coef_A_c_mean <- (smap_coef_A_c[1:31,]+smap_coef_A_c[33:63,]+smap_coef_A_c[65:95,])/3 

#Smap for Am

#Am E= 3
smap_res_Am_c <-  block_lnlp(block_data_c[,vars_c],
                           method = "s-map",
                           num_neighbors = 0, 
                           lib= segments_c, 
                           pred=segments_c,
                           theta = best_theta_c[2],
                           #columns = vars_c,
                           target_column = "Am",
                           silent = T,
                           save_smap_coefficients = T, 
                           first_column_time = F) 

#Time series of fluctuating interaction strength
smap_coef_Am_c <- as.data.frame(smap_res_Am_c$smap_coefficients[[1]])

smap_coef_Am_c_mean <- (smap_coef_Am_c[1:31,]+smap_coef_Am_c[33:63,]+smap_coef_Am_c[65:95,])/3 


#Smap for S
#S E= 3
#None causes S

smap_res_S_c <-  block_lnlp(block_data_c[,vars_c],
                             method = "s-map",
                             num_neighbors = 0, 
                             lib= segments_c, 
                             pred=segments_c,
                             theta = best_theta_c[3],
                             #columns = vars_c,
                             target_column = "S",
                             silent = T,
                             save_smap_coefficients = T, 
                             first_column_time = F) 

#Time series of fluctuating interaction strength
smap_coef_S_c <- as.data.frame(smap_res_S_c$smap_coefficients[[1]])

smap_coef_S_c_mean <- (smap_coef_S_c[1:31,]+smap_coef_S_c[33:63,]+smap_coef_S_c[65:95,])/3 

names(smap_coef_A_c_mean)   <- c("A", "Am", "S")
names(smap_coef_Am_c_mean)  <- c("A", "Am", "S")
names(smap_coef_S_c_mean)   <- c("A", "Am", "S")


all_coef_c <- rbind(
  melt(cbind(smap_coef_A_c_mean[,1:3], target=rep("A", nrow(smap_coef_A_c_mean)))),
  melt(cbind(smap_coef_Am_c_mean[,1:3], target=rep("Am", nrow(smap_coef_Am_c_mean)))),
  melt(cbind(smap_coef_S_c_mean[,1:3], target=rep("S", nrow(smap_coef_S_c_mean))))
)

palette3   <- c("darkred", "lightpink3",  "darkgreen")

ggplot(all_coef_c,aes(y=value, x=as.factor(as.character(target)), color=target))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_jitter(alpha=0.5)+
  geom_boxplot()+
  facet_wrap(.~variable, scales="free", drop = T)+
  scale_color_manual(values=palette3)+
  labs(x="Effect of", y="Interaction strength")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Plot 
all_coef_c$time <- rep(3:33, times=9)

ggplot(all_coef_c,aes(y=value, x=as.numeric(time), group=target, color=target))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line()+
  facet_wrap(.~variable, scales="free_y", drop = T, strip.position = "top")+
  scale_color_manual(values=palette3)+
  labs(x="Time", y="Interaction strength")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

aggdata2_c <- aggregate(all_coef_c$value, by = list(target=all_coef_c$target, 
                                                variable=all_coef_c$variable), FUN = mean, na.rm=T)

#Here is an overview of the positive-negative interactions
aggdata2_c$sign <- ifelse(aggdata2_c$x>0, "positive",
                        ifelse(aggdata2_c$x<0, "negative", "none"))


sign_c <- ggplot(aggdata2_c, aes(y = variable, x = target, fill = sign)) +
  geom_tile(color="white")+
  labs(color="Interaction", y="Target", x="Effect", title="Type of interactions")+
  scale_fill_manual(name="Interaction\ntype", values=c("#999999", "#E69F00"))+
  mytheme+
  theme(legend.position = "none")

freq_c<- ggplot(aggdata2_c, aes(x=abs(x),fill=sign)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  scale_color_manual(name= "Type", values=c("#999999", "#E69F00"))+
  scale_fill_manual(name= "Type", values=c("#999999", "#E69F00"))+
  labs(title="Skewness",x="Interaction Strenght", y = "Frequency")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = c(0.8, 0.8))


plot_grid(sign_c, freq_c)


#Dynamic stability 
#Create interaction matrix for each time point 
ac <- as.data.frame(cbind(time = seq(3,33), smap_coef_A_c_mean[,1:3], target=rep("A", nrow(smap_coef_A_c_mean))))
bc <- as.data.frame(cbind(time = seq(3,33), smap_coef_Am_c_mean[,1:3], target=rep("Am", nrow(smap_coef_Am_c_mean))))
cc <- as.data.frame(cbind(time = seq(3,33), smap_coef_S_c_mean[,1:3], target=rep("S", nrow(smap_coef_S_c_mean))))

#Interaction matrix
stab_data_c <- rbind.fill(ac,bc,cc) 
stab_data_c <- stab_data_c[,c("time","A", "Am","S","target")]

#NAs to 0 
stab_data_c[is.na(stab_data_c)] <- 1

#Order and split data into time points
stab_data_c_sort <- stab_data_c[order(stab_data_c$time),]
stab_c_split <- split(stab_data_c_sort[,2:4], stab_data_c_sort$time)


eigenvalues_c <- sapply(stab_c_split, function(x){
  max(Re(eigen(x)$values))
})

stability2_c <- melt(eigenvalues_c)
stability2_c$time <- rownames(stability2_c)

#Interaction indexes
time_mean_c <- sapply(stab_c_split, function(x){
  mean(abs(x[row(x)!=col(x)]))
})

time_median_c <- sapply(stab_c_split, function(x){
  median(abs(x[row(x)!=col(x)]))
})

time_max_c <- sapply(stab_c_split, function(x){
  max(abs(x[row(x)!=col(x)]))
})

weak_c <- as.vector(time_median_c)/as.vector(time_max_c) #Weak interactions index

mean.interactions_c <- melt(time_mean_c)
mean.interactions_c$time <- rownames(mean.interactions_c)

weak.interactions_c <- melt(weak_c)
weak.interactions_c$time <- rownames(stability2_c)

#Diversity index (exponential Shannon)
comp_split   <- split(competition, competition$replicate)
comp_shannon <- apply(competition[,3:5],1, function(x) {exp(diversity(x, index ="shannon"))})
comp_shannon_data <- cbind(competition[,1:2], comp_shannon)

library(dplyr)
grouped_shannon_c <- group_by(comp_shannon_data, day) %>%
  summarise(mean=mean(comp_shannon))

#Coefficient of variation 
cv <- function(x) {sd(x)/mean(x)}
grouped_abun_c  <- group_by(competition_ts, day, variable) %>%
  summarise(cv=cv(value+0.1)) %>%
  summarise(mean_cv=mean(cv))


#Plots 
stab_c <- ggplot(stability2_c,aes(y=value, x=as.numeric(time)))+
  geom_hline(yintercept = 1, linetype="dashed", color="darkred")+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y=expression(atop("Dynamic stability",paste(stable %<->% unstable))))+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

mean.int_c <- ggplot(mean.interactions_c,aes(y=value, x=as.numeric(time)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y="Mean interaction\nstrength")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

weak.int_c <- ggplot(weak.interactions_c,aes(y=value, x=as.numeric(time)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  geom_segment(mapping=aes(x=3, y=0.5, xend=3, yend=0.42), 
               size=0.7, color="black", 
               arrow = arrow(length = unit(0.4, "cm")))+
  annotate("text", label = "Dominance of\nweak interactions", x = 4.5, y = 0.45, 
           color = "black", angle=90)+
  labs(x="Time", y="Median:maximum\ninteraction strenght")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

comp_cv_c <- ggplot(grouped_abun_c,aes(y=mean_cv, x=as.numeric(day)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y="Coefficient of\nvariation")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

comp_diversity_c <- ggplot(grouped_shannon_c,aes(y=mean, x=as.numeric(day)))+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line()+
  labs(x="Time", y="Diversity index\nHill 1")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#See what causes dynamic stability

abun.morph_c <- as.data.frame(subset(dfc_comp, dfc_comp$variable=="Am"))

all_indexes_c <- as.data.frame(cbind(dyn.stab=stability2_c$value, mean.int=mean.interactions_c$value, 
                                   weak.int=weak.interactions_c$value, div=grouped_shannon_c$mean, 
                                   variation=grouped_abun_c$mean_cv, morph.abun=abun.morph_c$value))

all_indexes_norm_c <- as.data.frame(apply(all_indexes_c, 2, function(x) (x - min(x))/(max(x)-min(x))))

#Find best E
simplex_out <-  simplex(all_indexes_norm_c$dyn.stab, E=1:6, tau=1)

plot(simplex_out$E, simplex_out$rho, type="l", 
       xlab="Embedding Dimension (E)", ylab="Forecast Skill (rho)")

#E=1


#*********
###CCM###
#********

#CCM for tp=0 which is default
#Effect of mean interactions

ccm_output_stab_int_c <- ccm(all_indexes_norm_c, 
                           E=1,  
                           lib_column = "dyn.stab",
                           target_column = "mean.int",
                           random_libs = T,
                           num_samples = 10000,
                           silent=T)


ccm_output_stab_int_c$direction <- paste(ccm_output_stab_int_c$lib_column, "xmap to", ccm_output_stab_int_c$target_column)

ccm_ag_stab_int_c <- with(ccm_output_stab_int_c, aggregate(cbind(rho=rho),
                                                       list(direction=direction,
                                                            lib_column=lib_column,
                                                            lib_size=lib_size),
                                                       function(x) quantile(x, c(0.025,
                                                                                 pnorm(-1,0,1), 0.5, 
                                                                                 pnorm(1,0,1), 0.975), 
                                                                            na.rm=T)))

#Plot 
ccm_ag_stab_int_c$guide <- rep(0, length=nrow(ccm_ag_stab_int))
ccm_fig_stab_int_c <- ggplot(ccm_ag_stab_int_c, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Mean interaction strength")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_int_c


#Effect of weak interactions
ccm_output_stab_weak_c <- ccm(all_indexes_norm_c, 
                            E=1,  
                            lib_column = "dyn.stab",
                            target_column = "weak.int",
                            random_libs = T,
                            num_samples = 10000,
                            silent=T)


ccm_output_stab_weak_c$direction <- paste(ccm_output_stab_weak_c$lib_column, "xmap to", ccm_output_stab_weak_c$target_column)

ccm_ag_stab_weak_c <- with(ccm_output_stab_weak_c, aggregate(cbind(rho=rho),
                                                         list(direction=direction,
                                                              lib_column=lib_column,
                                                              lib_size=lib_size),
                                                         function(x) quantile(x, c(0.025,
                                                                                   pnorm(-1,0,1), 0.5, 
                                                                                   pnorm(1,0,1), 0.975), 
                                                                              na.rm=T)))

#Plot 
ccm_ag_stab_weak_c$guide <- rep(0, length=nrow(ccm_ag_stab_weak))
ccm_fig_stab_weak_c <- ggplot(ccm_ag_stab_weak_c, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Weak interactions")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_weak_c


#Plot 
ccm_ag_stab_div_c$guide <- rep(0, length=nrow(ccm_ag_stab_div))
ccm_fig_stab_div_c <- ggplot(ccm_ag_stab_div_c, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Diversity index")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_div_c


#Effect of coefficient of variation #no causation

#Effect of morphotype abundance
ccm_output_stab_morph_c <- ccm(all_indexes_norm_c, 
                             E=1,  
                             lib_column = "dyn.stab",
                             target_column = "morph.abun",
                             random_libs = T,
                             num_samples = 10000,
                             silent=T)


ccm_output_stab_morph_c$direction <- paste(ccm_output_stab_morph_c$lib_column, "xmap to", ccm_output_stab_morph_c$target_column)

ccm_ag_stab_morph_c<- with(ccm_output_stab_morph_c, aggregate(cbind(rho=rho),
                                                          list(direction=direction,
                                                               lib_column=lib_column,
                                                               lib_size=lib_size),
                                                          function(x) quantile(x, c(0.025,
                                                                                    pnorm(-1,0,1), 0.5, 
                                                                                    pnorm(1,0,1), 0.975), 
                                                                               na.rm=T)))

#Plot 
ccm_ag_stab_morph_c$guide <- rep(0, length=nrow(ccm_ag_stab_morph_c))
ccm_fig_stab_morph_c <- ggplot(ccm_ag_stab_morph_c, aes(x=lib_size, y=rho[,3]))+
  geom_line(aes(y=guide), color="grey", linetype="dashed")+
  geom_ribbon(aes(ymin=rho[,1], ymax=rho[,5], x=lib_size),alpha=0.1)+
  geom_line(color="darkred")+
  labs(y="Cross Map Skill (rho)", x="Library Size")+
  ggtitle("Morphotype abundance")+
  mytheme+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ccm_fig_stab_morph_c


plot_grid(ccm_fig_stab_int_c,ccm_fig_stab_weak_c, ccm_fig_stab_morph_c, ncol=3)


#Interaction strenghts 
test_theta_all <-
  block_lnlp(all_indexes_norm_c$dyn.stab,
             method = "s-map",
             num_neighbors = 0, #use any value < 1 for s-map
             theta = seq(1,4,0.1))

test_theta_all$theta[which.max(round(test_theta_all$rho,2))]

all_indexes_norm_c_selected <- all_indexes_norm_c[,c("dyn.stab", "mean.int", "weak.int", "morph.abun")]

smap_stab_c <-  block_lnlp(all_indexes_norm_c,
                         method = "s-map",
                         num_neighbors = 0,
                         theta = 1,
                         target_column = "dyn.stab",
                         silent = T,
                         save_smap_coefficients = T, 
                         first_column_time = F) 


#Time series of fluctuating interaction strength
smap_stab_dyn_c <- as.data.frame(smap_stab_c$smap_coefficients[[1]])

smap_stab_dyn_sub_c <- melt(smap_stab_dyn_c[,2:4])

contribution_c<- ggplot(smap_stab_dyn_sub_c, aes(x=variable, y=abs(value), fill=variable))+
  geom_boxplot()+
  ylab("Influences on stability")+
  scale_x_discrete(labels = c("Mean interaction\ninteraction strenght", 'Dominance of\nweak interactions', "Morphotype\nabundance"))+
  scale_fill_manual(values = c("#999999", "#E69F00", "#999999"))+
  mytheme+
  theme(legend.position= "none",
        axis.title.x = element_blank(), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))


#Combined Interaction plot 
plot_grid(stab_c, mean.int_c, weak.int_c, comp_cv_c, comp_diversity_c, contribution_c,
          ncol=2)


#Draw the networks

#Create incidence matrices 
#signif
#signif_tp
incidences2_c <- (signif_c %>% separate(direction, 
                                    c("variable", "rest","rest", "target")))[, c(1,4,6)]

#add the positive-negative information 
network.data2_c <- (merge(incidences2_c, aggdata2_c, by=c("variable", "target"), all=T))
#names(network.data2_c) <- c("to", "from", "significance", "type")

palette6  <- c("grey25", "lightpink3", "skyblue3", "darkgreen", "midnightblue", "darkred")
#Subset the significant causations for plotting purposes 
#network.data.subset2_c <- subset(network.data2_c, network.data2_c$significance==1)


#Subset the significant causations for plotting purposes 
#network.data.subset2_c <- subset(network.data2_c, network.data2_c$significance==1)

network.data2_c$variable <-mapvalues(network.data2_c$variable, from = c("A", "Am","S"), 
                                  to=c("3A","2Am","1S"))
network.data2_c$target <-mapvalues(network.data2_c$target, from = c("A", "Am","S"), 
                                   to=c("3A", "2Am","1S" ))
palette6   <- c("darkgreen","lightpink3","darkred")


net_c <- ggplot(data = network.data2_c, aes(from_id = target, to_id = variable)) +
  geom_net(layout.alg = 'circle', labelon = TRUE, 
           size = 11, directed = TRUE, vjust = 0.5, 
           arrowsize = 1, arrowgap = 0.07,
           selfloops = F, singletons= TRUE, labelcolour = "white", 
           ealpha = .7, fontsize = 5, curvature=0.1,
           aes(colour = target,linewidth = rescale(abs(x), c(0.4,6))))+ 
  theme_net()+
  scale_colour_manual(values=palette6)+
  scale_linetype_manual(values=c("twodash", "solid"))+
  theme(panel.background = element_rect(color = "grey"))+
  theme(legend.position = "none")

#Combined plots

plot_grid(net, net_c, labels = c("A", "B"), 
          align = 'h')

#Test
palette4 <- c("darkred", "midnightblue", "darkgreen" , "grey25")
pairwise.data <- read.table("competitive_effect_exp.csv", header=T,sep=";",dec=",")

pairwise <- ggplot(data = pairwise.data, aes(x=EffectOf, y=Effect, color=EffectOf)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot()+
  labs(x="Effect of", y="Interaction strenght")+
  facet_wrap(~Variable, drop=T, ncol=2, scales="free")+
  scale_color_manual(values=palette4)+
  mytheme+
theme(legend.position= "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

#Combined plots with without predator 
plot_grid(stab, stab_c, labels=c("Predation", "Competition"))
plot_grid(mean.int, mean.int_c, labels=c("Predation", "Competition"))
plot_grid(weak.int, weak.int_c, labels=c("Predation", "Competition"))
plot_grid(pred_diversity, comp_diversity_c, labels=c("Predation", "Competition"))
plot_grid(contribution, contribution_c, labels=c("Predation", "Competition"))
plot_grid(sign, sign_c, labels=c("Predation", "Competition"))
plot_grid(pred_cv, comp_cv_c, labels=c("Predation", "Competition"))


stability$var <- rep("predation", times=nrow(stability))
stability2_c$var <- rep("competition", times=nrow(stability))

stability_new <-as.data.frame(rbind(stability, stability2_c))

stab_c <- ggplot(stability_new,aes(y=value, x=as.numeric(time), group=var))+
  geom_hline(yintercept = 1, linetype="dashed", color="red", size=1)+
  geom_vline(xintercept = 14, linetype="dashed", color="grey")+
  geom_vline(xintercept = 25, linetype="dashed", color="grey")+
  geom_line(size=1, aes(linetype=var))+
  labs(x="Time", y="Stability")+
  mytheme+
  theme(legend.position="none", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
