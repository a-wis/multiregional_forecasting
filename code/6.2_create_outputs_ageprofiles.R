#############################
## Plotting age profiles
#############################
#creating plots and tables for the manuscript
#A unified framework for probabilistic forecasting subnational populations
#(C) anonymised (2020) 

# the code first transforms the results from stan and data into handy data frames, then produces required summaries, then plots using ggplot2+ggfan+ggpubr. Saving of the plots has been disabled.

# reading libraries #### 
require(dplyr)
require(tidyr)
require(forcats)
require(readr)
require(rstan)
require(ggplot2)
require(ggfan)
require(ggpubr)
require(stats)
require(ggh4x)
#the code also uses mapvalues funcation from plyr package

#od
aodd = input_au %>% group_by(Reg_O,Reg_D,year,S,A) %>% #imi for regions
  summarise(dat_x=(Flow/x)) %>%
  ungroup() %>%
  mutate(Reg_O=as.factor(Reg_O), Reg_D=as.factor(Reg_D)) %>%
  mutate(Reg_O=plyr::mapvalues(Reg_O,1:8,Regn),
         Reg_D=plyr::mapvalues(Reg_D,1:8,Regn)) %>% 
  rename(Sex=S, Age=A) 

aodestl = t_od %>% mutate(Reg_O=as.factor(Reg_O), 
                           Reg_D=as.factor(Reg_D)) %>%
  group_by(Reg_D,Reg_O,year,S,A) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>% 
  mutate(Reg_O=plyr::mapvalues(Reg_O,1:8,Regn),Reg_D=plyr::mapvalues(Reg_D,1:8,Regn)) %>%
  rename(Sex=S, Age=A) %>%
  left_join(aodd) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile), 
         Age=plyr::mapvalues(Age,levels(Age),agelab[-18]),
         x=ifelse(Sex=="Males",-x,x),
         dat_x=ifelse(Sex=="Males",-dat_x,dat_x)) 
aodestl$dat_x[aodestl$year==2021]<-aodestl$dat_x[aodestl$year==2011]
aodestl$dat_x[aodestl$year==2016]<-aodestl$dat_x[aodestl$year==2011]

# internal migration age####
plot.od.age<-function(ye){
  p=ggplot(data=aodestl %>% filter(year==ye-5), #lex0.04 %>% select(-A)
           aes(x=Age, y=x, group=Sex, color=Sex, quantile=quantile)) +
    geom_fan(aes(group=Sex), intervals=c(0.5,0.8,0.95)) +
    # geom_col(aes(x=Age, y=dat_x, color=Sex, quantile=NULL), width=1, alpha=0) +
    geom_line(aes(x=Age, y=dat_x, linetype=Sex, quantile=NULL), size=1.02,alpha=1,colour="black") + #colour="black"
    scale_fill_gradient(name="Estimates\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + 
    scale_linetype(name="Sex", label=c(paste0("Males ",if (ye>2016) 2016 else ye),paste0("Females ",if (ye>2016) 2016 else ye))) +
    scale_y_continuous(labels=function(x) abs(x), name="Out-migration probability") +
    scale_x_discrete(expand=expand_scale(mult = c(0, 0))) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0)))) +
    facet_grid( Reg_O ~ Reg_D, scales = "free_x", switch="y") + theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), 
          axis.text.x  = element_text(angle=90, vjust=0.5), 
          axis.text.y  = element_text(size=6), 
          legend.position = "bottom") +
    coord_flip()
  ggsave(paste0("plots/od_age",ye,".pdf"),plot=p,device="pdf",width = 4, height = 4, scale=2.5)
}
plot.od.age(2016)
plot.od.age(2026)

#fertility
afert16 = t_f16 %>% group_by(Reg,A,year) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>% 
  left_join(births_8116) %>%
  select(-age,-state_name) %>%
  mutate(Reg=plyr::mapvalues(Reg,1:8,Regn)) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile)) %>%
  rename(Age=A)
for (i in 2017:2026){
  afert16$ERP[afert16$year==i]<-afert16$ERP[afert16$year==2016]
  afert16$births[afert16$year==i]<-afert16$births[afert16$year==2016]    
}


# fertility age ####
plot.tfr.age<-function(ye){
  p=ggplot(data=afert16 %>% filter(year==ye), 
           aes(x=Age, y=x, group=1, quantile=quantile)) +
    geom_fan(intervals=c(0.5,0.8,0.95)) +
    # geom_col(aes(x=Age, y=dat_x, color=Sex, quantile=NULL), width=1, alpha=0) +
    geom_line(aes(x=Age, y=births/ERP, quantile=NULL, linetype="Data"),  size=1.02,alpha=1,colour="black") + #colour="black"
    geom_point(aes(x=Age, y=births/ERP, quantile=NULL, shape="Data"), alpha=1, size=2,colour="black") +
    scale_fill_gradient(name="Estimates\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + 
    scale_y_continuous(labels=function(x) abs(x), name="Fertility rate") +
    scale_x_discrete(expand=expand_scale(mult = c(.02, 0.02))) +
    scale_linetype(name="", label=paste0("Data ",if (ye>2016) 2016 else ye)) +
    scale_shape(name="", label=paste0("Data ",if (ye>2016) 2016 else ye)) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
           linetype=guide_legend(override.aes = list(linetype=1, shape=16))) +
    facet_wrap( Reg ~ ., scales = "fixed", ncol=4) + theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), axis.text.x  = element_text(angle=90, vjust=0.5), legend.position = "bottom") 
  ggsave(paste0("./Plots/tfr_age",ye,".pdf"),plot=p,device="pdf",width = 5, height = 3, scale=2.2)
}

#plot.tfr.age(2016)

#immigration age ####
aimi44l = t_imi44 %>% group_by(Reg,year,S,A) %>%
  summarise(`0.025`=quantile(sqrt(x),c(.025)),
            `0.1`=quantile(sqrt(x),c(.1)),
            `0.25`=quantile(sqrt(x),c(.25)),
            `0.5`=quantile(sqrt(x),c(.5)),
            `0.75`=quantile(sqrt(x),c(.75)),
            `0.9`=quantile(sqrt(x),c(.9)),
            `0.975`=quantile(sqrt(x),c(.975))) %>%
  ungroup() %>%
  left_join(intmigr_8116) %>%
  mutate(Reg=plyr::mapvalues(Reg,1:8,Regn)) %>% rename(Sex=S, Age=A) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile)) 
prepare.imi.plot<-function(inp){
  inp$x[inp$Sex=="M"]<- -inp$x[inp$Sex=="M"]
  inp$emig[inp$Sex=="M"]<- -inp$emig[inp$Sex=="M"]
  inp$immig[inp$Sex=="M"]<- -inp$immig[inp$Sex=="M"]
  for (i in 2017:2026){
    inp$ERP[inp$year==i]<-inp$ERP[inp$year==2016]
    inp$emig[inp$year==i]<-inp$emig[inp$year==2016]
    inp$immig[inp$year==i]<-inp$immig[inp$year==2016]    
  }
  return(inp)
}
aimi44l=prepare.imi.plot(aimi44l)

plot.imi.age<-function(inp=aimi02l,ye,nam){
  p=ggplot(data=inp %>% filter(year==ye), #lex0.04 %>% select(-A)
           aes(x=Age, y=x, group=Sex, color=Sex, quantile=quantile)) +
    geom_fan(aes(group=Sex), intervals=c(0.5,0.8,0.95)) +
    # geom_col(aes(x=Age, y=dat_x, color=Sex, quantile=NULL), width=1, alpha=0) +
    geom_line(aes(x=Age, y=sign(immig)*sqrt(abs(immig)), linetype=Sex, quantile=NULL), size=1.02,alpha=1,colour="black") + #colour="black"
    scale_fill_gradient(name="Estimates\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + 
    scale_linetype(name="Sex", label=c(paste0("Males ",if (ye>2016) 2016 else ye),paste0("Females ",if (ye>2016) 2016 else ye))) +
    scale_y_continuous(labels=function(x) abs(x), 
                       name="Immigration counts (square root)", 
                       limits = c(-200,200)) +
    scale_x_discrete(expand=expand_scale(mult = c(0, 0))) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0)))) +
    facet_wrap( Reg ~ ., scales = "fixed", ncol=4) + theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), 
          axis.text.x  = element_text(angle=90, vjust=0.5), 
          legend.position = "bottom", 
          legend.text=element_text(size=12)) +
    coord_flip()
  ggsave(paste0("plots/imi_",nam,"_age",ye,".pdf"),plot=p,device="pdf",width = 5, height = 3, scale=2.2)
}
plot.imi.age(aimi44l,2016,"031")

# emigration ####
prepare.emi.plot<-function(inp=t_emi02){
  ret= inp %>% group_by(Reg,year,S,A) %>%
    summarise(`0.025`=quantile(sqrt(x),c(.025)),
              `0.1`=quantile(sqrt(x),c(.1)),
              `0.25`=quantile(sqrt(x),c(.25)),
              `0.5`=quantile(sqrt(x),c(.5)),
              `0.75`=quantile(sqrt(x),c(.75)),
              `0.9`=quantile(sqrt(x),c(.9)),
              `0.975`=quantile(sqrt(x),c(.975))) %>%
    ungroup() %>%
    left_join(intmigr_8116) %>%
    mutate(Reg=plyr::mapvalues(Reg,1:8,Regn)) %>% rename(Sex=S, Age=A) %>%
    gather("quantile","x",`0.025`:`0.975`) %>%
    mutate(quantile=as.numeric(quantile)) 
  ret$x[ret$Sex=="M"]<- -ret$x[ret$Sex=="M"]
  ret$emig[ret$Sex=="M"]<- -ret$emig[ret$Sex=="M"]
  ret$immig[ret$Sex=="M"]<- -ret$immig[ret$Sex=="M"]
  for (i in 2017:2026){
    ret$ERP[ret$year==i]<-ret$ERP[ret$year==2016]
    ret$emig[ret$year==i]<-ret$emig[ret$year==2016]
    ret$immig[ret$year==i]<-ret$immig[ret$year==2016]    
  }
  return(ret)
}

aemi02l = prepare.emi.plot(t_emi44)
# aemi44l = prepare.emi.plot(t_emi44)


plot.emi.age<-function(inp,ye,nam){
  p=ggplot(data=inp %>% filter(year==ye), #lex0.04 %>% select(-A)
           aes(x=Age, y=x, group=Sex, color=Sex, quantile=quantile)) +
    geom_fan(aes(group=Sex), intervals=c(0.5,0.8,0.95)) +
    # geom_col(aes(x=Age, y=dat_x, color=Sex, quantile=NULL), width=1, alpha=0) +
    geom_line(aes(x=Age, y=sign(emig)*sqrt(abs(emig)/abs(ERP)), linetype=Sex, quantile=NULL), size=1.02,alpha=1,colour="black") + #colour="black"
    scale_fill_gradient(name="Estimates\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + 
    scale_linetype(name="Data", label=c(paste0("Females ",if (ye>2016) 2016 else ye),paste0("Males ",if (ye>2016) 2016 else ye))) +
    scale_y_continuous(labels=function(x) abs(x), 
                       name="Emigration rate (square root)",
                       limits=c(-.4,0.4)) +
    scale_x_discrete(expand=expand_scale(mult = c(0, 0))) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0)))) +
    facet_wrap( Reg ~ ., scales = "fixed", ncol=4) + theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), axis.text.x  = element_text(angle=90, vjust=0.5), legend.position = "bottom", legend.text=element_text(size=12)) +
    coord_flip()
  ggsave(paste0("plots/emi_",nam,"_age",ye,".pdf"),plot=p,device="pdf",width = 5, height = 3, scale=2.2)
}
# for (i in c(2006,2011,2016,2026)) plot.imi.age()
for (i in c(2006,2016,2026)) plot.emi.age(aemi02l,i,"02")



####
#mortality ####
####
prepare.mort.plot<-function(inp=t_mr14){
  ret= inp %>% group_by(Reg,year,S,A) %>%
    summarise(`0.025`=quantile(log10(x),c(.025)),
              `0.1`=quantile(log10(x),c(.1)),
              `0.25`=quantile(log10(x),c(.25)),
              `0.5`=quantile(log10(x),c(.5)),
              `0.75`=quantile(log10(x),c(.75)),
              `0.9`=quantile(log10(x),c(.9)),
              `0.975`=quantile(log10(x),c(.975))) %>%
    ungroup() %>%
    left_join(deaths9311) %>%
    mutate(Reg=plyr::mapvalues(Reg,1:8,Regn)) %>% rename(Sex=S, Age=A, dat_x=x) %>%
    gather("quantile","x",`0.025`:`0.975`) %>%
    mutate(quantile=as.numeric(quantile)) 
  ret$x[ret$Sex=="F"]<- -ret$x[ret$Sex=="F"]
  ret$dat_x[ret$Sex=="F"]<- -ret$dat_x[ret$Sex=="F"]
  #take log10 of data
  
  for (i in 2017:2026){
    ret$ERP[ret$year==i]<-ret$ERP[ret$year==2011]
    ret$dat_x[ret$year==i]<-ret$dat_x[ret$year==2011]
  }
  return(ret)
}

tmort14=prepare.mort.plot(t_mr14)

plot.mort.age<-function(inp,ye,nam){
  p=ggplot(data=inp %>% filter(year==ye), #lex0.04 %>% select(-A)
           aes(x=Age, y=x, group=Sex, color=Sex, quantile=quantile)) +
    geom_fan(aes(group=Sex), intervals=c(0.5,0.8,0.95)) +
    # geom_col(aes(x=Age, y=dat_x, color=Sex, quantile=NULL), width=1, alpha=0) +
    geom_line(aes(x=Age, y=sign(dat_x)*log10(abs(dat_x)/abs(ERP)), linetype=Sex, quantile=NULL), size=0.8,alpha=1,colour="black") + #colour="black"
    scale_fill_gradient(name="Estimates\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + 
    scale_linetype(name="Data", label=c(paste0("Females ",if (ye>2011) 2011 else ye),paste0("Males ",if (ye>2016) 2016 else ye))) +
    scale_y_continuous(labels=function(x) -abs(x), 
                       name="Mortality rate (common logarithm)") +
    scale_x_discrete(expand=expand_scale(mult = c(0, 0))) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0)))) +
    facet_wrap( Reg ~ ., scales = "fixed", ncol=4) + theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), axis.text.x  = element_text(angle=90, vjust=0.5), legend.position = "bottom", legend.text=element_text(size=12)) +
    coord_flip()
  ggsave(paste0("plots/mortL_",nam,"_age",ye,".pdf"),plot=p,device="pdf",width = 5, height = 3, scale=2.2)
}

for (i in c(2006,2011,2016,2026)) plot.mort.age(tmort14,i,"02")
