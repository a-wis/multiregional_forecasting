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

source("code/1_data_processing.R")
source("code/5.0_functions_forecasts_transformations.R")


#transformations #####
# mat_m=gen.forecast(inp = fit_m0124,S=2,R=8,Fh = 15,Ti=19,Ag=18,fit.too =T)
# mat_f=gen.forecast(inp = fit_f0016,S=1,R=8,Fh = 15,Ti=31,Ag=7,fit.too = T)
# mat_im=gen.forecast(inp = fit_imi0044,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =T)
# mat_em=gen.forecast(inp = fit_emi0044,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =T)
# mat_od=gen.forecast(inp = fit05_01,S=2,R=56,Fh = 3,Ti=6,Ag=17,fit.too =T)
mat_m=gen.forecast(inp = fit_m02,S=2,R=8,Fh = 15,Ti=19,Ag=18,fit.too =T)
mat_f=gen.forecast(inp = fit_f01,S=1,R=8,Fh = 15,Ti=31,Ag=7,fit.too = T)
mat_im=gen.forecast(inp = fit_imi031,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =T)
mat_em=gen.forecast(inp = fit_emi03,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =T)
mat_od=gen.forecast(inp = fit08_04,S=2,R=56,Fh = 3,Ti=6,Ag=17,fit.too =T)

#creating population projection (or has to be read in, if previously saved)
source("code/create_projection.R")

## mortality ####
t_m14=transform.mort.lex(mat_m)
t_mr14=transform.df(mat_obj = mat_m, dat = deaths9311,yd0=1993,ydl=2011,yf0=2012,yfl=2026,pop.comp = "M",interval = 1)

###transforming mortality data
lexd=life.tab.d()
lexd=apply(lexd,c(3:5),diag)
lexd=as.data.frame.table(lexd)
lexd = lexd %>% rename(Reg=Var1, A=Var2, S=Var3, year=Var4,x=Freq) %>%
  mutate(Reg=plyr::mapvalues(Reg,levels(Reg),Regn),
         A=plyr::mapvalues(A,levels(A),levels(deaths9311$A)),
         S=plyr::mapvalues(S,levels(S),levels(deaths9311$S)),
         year=as.integer(1992+as.integer(plyr::mapvalues(year,levels(year),c(1993:2011)))))

## transforming rest of the results into df ####
t_f16 = transform.df(mat_obj = mat_f, dat = births_8116,pop.comp = "F")
t_imi44 = transform.df(mat_obj = mat_im, dat = intmigr_8116,pop.comp = "I")
t_emi44 = transform.df(mat_obj = mat_em, dat = intmigr_8116,pop.comp = "E")
t_od = transform.df(mat_obj = mat_od, yd0=1981, ydl=2006, yf0=2007, yfl=2021, dat = input_au,pop.comp = "OD")



#preparation of data for plotting ####
#Preparation: internal migration ####
##by sex
todd = input_au %>% group_by(Reg_O,Reg_D,Corridor,year,S) %>% #imi for regions
  summarise(dat_x=log10(mean(Flow/x))) %>%
  ungroup() %>%
  mutate(Reg_O=as.factor(Reg_O), Reg_D=as.factor(Reg_D)) %>%
  mutate(Reg_O=plyr::mapvalues(Reg_O,1:8,Regn),Reg_D=plyr::mapvalues(Reg_D,1:8,Regn)) %>% rename(Sex=S)

tod3111l = t_od %>% mutate(Reg_O=as.factor(Reg_O), Reg_D=as.factor(Reg_D)) %>%
  group_by(Reg_O,Reg_D,year,S,Iter) %>%
  summarise(x=log10(mean(x))) %>% 
  ungroup() %>%
  group_by(Reg_O,Reg_D,year,S) %>%
  # summarise(q25=quantile(x,c(.1)),q50=quantile(x,c(.5)),q75=quantile(x,c(.9))) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>% 
  mutate(Reg_O=plyr::mapvalues(Reg_O,1:8,Regn),Reg_D=plyr::mapvalues(Reg_D,1:8,Regn)) %>%
  rename(Sex=S) %>%
  left_join(todd) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile),
         levels=10^(x)) 

##total
todd = input_au %>% group_by(Reg_O,Reg_D,Corridor,year) %>% #S removed
  summarise(dat_x=log10(mean(Flow/x))) %>%
  ungroup() %>%
  mutate(Reg_O=as.factor(Reg_O), Reg_D=as.factor(Reg_D)) %>%
  mutate(Reg_O=plyr::mapvalues(Reg_O,1:8,Regn),Reg_D=plyr::mapvalues(Reg_D,1:8,Regn))

tod3111l = t_od %>% mutate(Reg_O=as.factor(Reg_O), Reg_D=as.factor(Reg_D)) %>%
  group_by(Reg_O,Reg_D,year,Iter) %>%
  summarise(x=log10(mean(x))) %>% 
  ungroup() %>%
  group_by(Reg_O,Reg_D,year) %>%
  # summarise(q25=quantile(x,c(.1)),q50=quantile(x,c(.5)),q75=quantile(x,c(.9))) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>% 
  mutate(Reg_O=plyr::mapvalues(Reg_O,1:8,Regn),Reg_D=plyr::mapvalues(Reg_D,1:8,Regn)) %>%
  left_join(todd) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile),
         levels=10^(x)) 

# levels (NOT USED)
tod3111 = t_od %>% mutate(Reg_O=as.factor(Reg_O), Reg_D=as.factor(Reg_D)) %>%
  group_by(Reg_O,Reg_D,year,S,Iter) %>%
  summarise(x=(mean(x))) %>% 
  ungroup() %>%
  group_by(Reg_O,Reg_D,year,S) %>%
  # summarise(q25=quantile(x,c(.1)),q50=quantile(x,c(.5)),q75=quantile(x,c(.9))) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>% 
  mutate(Reg_O=plyr::mapvalues(Reg_O,1:8,Regn),Reg_D=plyr::mapvalues(Reg_D,1:8,Regn)) %>%
  rename(Sex=S) %>%
  left_join(todd) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile),
         levels=exp(x)) 




#Preparation: fertility ####
#TFR total for all regions by year 
tfrdt = births_8116 %>% group_by(age,year) %>% #total TFR for AU
  summarise(births=sum(births),ERP=sum(ERP)) %>% 
  ungroup() %>%
  group_by(year) %>%
  summarise(x=sum(births/ERP)*5) %>% 
  ungroup()
#fertility by age-region-year
# temp16av=temp16 %>% group_by(Reg,year,A) %>%
#   summarise(q25=quantile(x,c(.25)),q50=quantile(x,c(.5)),q75=quantile(x,c(.75))) %>%
#   ungroup()

#TFR - data by region-year
tfrd = births_8116 %>% group_by(Reg,year) %>% #TFR for regions
  summarise(dat_x=sum(births/ERP)*5) %>%
  ungroup() 
#TFR - summarising results and joining with data 
tfr16 = t_f16 %>% group_by(Reg,year,Iter) %>%
  summarise(x=sum(x)*5) %>% 
  ungroup() %>%
  group_by(Reg,year) %>%
  # summarise(q25=quantile(x,c(.1)),q50=quantile(x,c(.5)),q75=quantile(x,c(.9))) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>% 
  left_join(tfrd) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile), Reg=plyr::mapvalues(Reg,1:8,Regn))


#Preparation: mortality ####
#lex
lexd0 = lexd %>% group_by(Reg,A,S,year) %>% #imi for regions
  filter(A=="0-4") %>% 
  rename(dat_x=x) %>%
  ungroup()
lex0.14 = t_m14 %>% filter(A=="0-4") %>% 
  group_by(Reg,A,S,year) %>%
  # summarise(q25=quantile(x,c(.1)),q50=quantile(x,c(.5)),q75=quantile(x,c(.9))) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>%
  left_join(lexd0) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile),
         S = case_when(S=="F" ~ "Female", S=="M" ~ "Male")) %>%
  rename(Sex=S)



#Preparation: immigration ####
trans.imig.tot<-function(inp=t_imi04, FUN=log){
  timil = intmigr_8116 %>% group_by(Reg,year,S) %>% #imi for regions
    summarise(dat_x=FUN(sum(immig))) %>%
    ungroup()
  
  out = inp %>% group_by(Reg,year,S,Iter) %>%
    summarise(x=FUN(sum(x))) %>% 
    ungroup() %>%
    group_by(Reg,year,S) %>%
    summarise(`0.025`=quantile(x,c(.025)),
              `0.1`=quantile(x,c(.1)),
              `0.25`=quantile(x,c(.25)),
              `0.5`=quantile(x,c(.5)),
              `0.75`=quantile(x,c(.75)),
              `0.9`=quantile(x,c(.9)),
              `0.975`=quantile(x,c(.975))) %>%
    ungroup() %>%
    left_join(timil) %>%
    gather("quantile","x",`0.025`:`0.975`) %>%
    mutate(quantile=as.numeric(quantile)) %>%
    mutate(Reg=plyr::mapvalues(Reg,1:8,Regn), S = case_when(S=="F" ~ "Female", S=="M" ~ "Male")) %>% rename(Sex=S)
  return(out)
}
timi031l=trans.imig.tot(t_imi44,FUN=log10)
# timi44l10=trans.imig.tot(t_imi44,FUN=log10)

t_imid = intmigr_8116 %>% group_by(year) %>% #imi for regions
  summarise(dat_x=sum(immig)) %>%
  ungroup()
t_imitot=t_imi44 %>% group_by(year,Iter) %>%
  summarise(x=sum(x)) %>% 
  ungroup() %>%
  group_by(year) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>%
  left_join(t_imid) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile)) 


#Preparation: emigration ####
trans.emig.tot<-function(inp,FUN=log10){
  temil = intmigr_8116 %>% group_by(Reg,year,S) %>% #imi for regions
    summarise(dat_x=FUN(sum(emig/ERP))) %>%
    ungroup()
  out = inp %>% group_by(Reg,year,S,Iter) %>%
    summarise(x=FUN(sum(x))) %>% 
    ungroup() %>%
    group_by(Reg,year,S) %>%
    summarise(`0.025`=quantile(x,c(.025)),
              `0.1`=quantile(x,c(.1)),
              `0.25`=quantile(x,c(.25)),
              `0.5`=quantile(x,c(.5)),
              `0.75`=quantile(x,c(.75)),
              `0.9`=quantile(x,c(.9)),
              `0.975`=quantile(x,c(.975))) %>%
    ungroup() %>%
    left_join(temil) %>%
    gather("quantile","x",`0.025`:`0.975`) %>%
    mutate(quantile=as.numeric(quantile)) %>%
    mutate(Reg=plyr::mapvalues(Reg,1:8,Regn), S = case_when(S=="F" ~ "Female", S=="M" ~ "Male")) %>% rename(Sex=S)
  return(out)
}
temi02l=trans.emig.tot(t_emi44,FUN=log10)
# temi44l10=trans.emig.tot(t_emi44,FUN=log10)

#joining immigration and emigration into one data.f
migrl10 = bind_rows("Immigration counts"=timi031l,"Emigration rates"=temi02l,.id = "component")



#Tranformation: population forecasts####
#Preparation: population age/time/region####
Pf.df=as.data.frame.table(Pf_test)
Pf.df1 = Pf.df  %>% rename(Iter=Var1, Reg=Var2, A=Var3, S=Var4, year=Var5,x=Freq) %>%
  mutate(Iter=as.character(plyr::mapvalues(Iter,levels(Iter),1:1000)),
         Reg=plyr::mapvalues(Reg,levels(Reg),1:8),
         A=plyr::mapvalues(A,levels(A),levels(ERP_pop$A)),
         S=plyr::mapvalues(S,levels(S),levels(ERP_pop$S)),
         year=as.integer(2011+5*as.integer(plyr::mapvalues(year,levels(year),seq(2016,2026,5))))) %>%
  group_by(Reg,year,S,A) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  full_join(ERP_pop) %>%
  mutate(quantile=as.numeric(quantile),Reg=plyr::mapvalues(Reg,1:8,Regn), Sex=plyr::mapvalues(S,levels(S),c("Females","Males"))) %>% rename(Age=A) 
Pf.df1$ERP[Pf.df1$year==2026]<-Pf.df1$ERP[Pf.df1$year==2016]
Pf.df1$ERP[Pf.df1$year==2021]<-Pf.df1$ERP[Pf.df1$year==2016]
Pf.df1$x[Pf.df1$Sex=="Males"]<- -Pf.df1$x[Pf.df1$Sex=="Males"]

#manipulating data
ERP_pop0 = ERP_pop %>% mutate(Reg=plyr::mapvalues(Reg,1:8,Regn), Sex=plyr::mapvalues(S,levels(S),c("Females","Males"))) %>% rename(Age=A) 
ERP_pop0$ERP[ERP_pop0$Sex=="Males"]<- -ERP_pop0$ERP[ERP_pop0$Sex=="Males"]

ERP_pop2 = ERP_pop %>% group_by(Reg,year) %>% #imi for regions
  summarise(x=sum(ERP)) %>%
  ungroup()

Pf.df2 = Pf.df  %>% rename(Iter=Var1, Reg=Var2, A=Var3, S=Var4, year=Var5,x=Freq) %>%
  mutate(Iter=as.character(plyr::mapvalues(Iter,levels(Iter),1:1000)),
         Reg=plyr::mapvalues(Reg,levels(Reg),1:8),
         A=plyr::mapvalues(A,levels(A),levels(ERP_pop$A)),
         S=plyr::mapvalues(S,levels(S),levels(ERP_pop$S)),
         year=as.integer(2011+5*as.integer(plyr::mapvalues(year,levels(year),seq(2016,2026,5))))) %>%
  group_by(Iter,Reg,year) %>%
  summarise(x=sum(x)) %>%
  ungroup() %>%
  group_by(Reg,year) %>%
  summarise(q25=quantile(x,c(.1)),q50=quantile(x,c(.5)),q75=quantile(x,c(.9))) %>%
  ungroup() %>%
  full_join(ERP_pop2) %>%
  mutate(Reg=plyr::mapvalues(Reg,1:8,Regn)) 

ERP_pop3 = ERP_pop %>% group_by(Reg,year,S) %>% #imi for regions
  summarise(ERP=sum(ERP)) %>%
  ungroup()
ERP_pop30 = ERP_pop3 %>% mutate(Reg=plyr::mapvalues(Reg,1:8,Regn), Sex=plyr::mapvalues(S,levels(S),c("Females","Males"))) 

#Preparation: population region/time####
Pf.df3 = Pf.df  %>% rename(Iter=Var1, Reg=Var2, A=Var3, S=Var4, year=Var5,x=Freq) %>%
  mutate(Iter=as.character(plyr::mapvalues(Iter,levels(Iter),1:1000)),
         Reg=plyr::mapvalues(Reg,levels(Reg),1:8),
         A=plyr::mapvalues(A,levels(A),levels(ERP_pop$A)),
         S=plyr::mapvalues(S,levels(S),levels(ERP_pop$S)),
         year=as.integer(2011+5*as.integer(plyr::mapvalues(year,levels(year),seq(2016,2026,5))))) %>%
  group_by(Iter,Reg,S,year) %>%
  summarise(x=sum(x)) %>%
  ungroup() %>%
  group_by(Reg,S,year) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>%
  full_join(ERP_pop3) %>%
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile),Reg=plyr::mapvalues(Reg,1:8,Regn), Sex=plyr::mapvalues(S,levels(S),c("Females","Males"))) 

Pf.df3$quantile[Pf.df3$year==2011]<-Pf.df3$quantile[Pf.df3$year==2016]
Pf.df3$x[Pf.df3$year==2011]<-Pf.df3$ERP[Pf.df3$year==2011]

#Preparation: population total over time ####
Pf.df4 = Pf.df  %>% rename(Iter=Var1, Reg=Var2, A=Var3, S=Var4, year=Var5,x=Freq) %>%
  mutate(Iter=as.character(plyr::mapvalues(Iter,levels(Iter),1:1000)),
         Reg=plyr::mapvalues(Reg,levels(Reg),1:8),
         A=plyr::mapvalues(A,levels(A),levels(ERP_pop$A)),
         S=plyr::mapvalues(S,levels(S),levels(ERP_pop$S)),
         year=as.integer(2011+5*as.integer(plyr::mapvalues(year,levels(year),seq(2016,2026,5))))) %>%
  group_by(Iter,year) %>%
  summarise(x=sum(x)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(`0.025`=quantile(x,c(.025)),
            `0.1`=quantile(x,c(.1)),
            `0.25`=quantile(x,c(.25)),
            `0.5`=quantile(x,c(.5)),
            `0.75`=quantile(x,c(.75)),
            `0.9`=quantile(x,c(.9)),
            `0.975`=quantile(x,c(.975))) %>%
  ungroup() %>% 
  gather("quantile","x",`0.025`:`0.975`) %>%
  mutate(quantile=as.numeric(quantile))
ERP_pop4 = ERP_pop %>% group_by(year) %>% #imi for regions
  summarise(ERP=sum(ERP)) %>%
  ungroup() 

Pf.df40 = Pf.df  %>% rename(Iter=Var1, Reg=Var2, A=Var3, S=Var4, year=Var5,x=Freq) %>%
  mutate(Iter=as.character(plyr::mapvalues(Iter,levels(Iter),1:1000)),
         Reg=plyr::mapvalues(Reg,levels(Reg),1:8),
         A=plyr::mapvalues(A,levels(A),levels(ERP_pop$A)),
         S=plyr::mapvalues(S,levels(S),levels(ERP_pop$S)),
         year=as.integer(2011+5*as.integer(plyr::mapvalues(year,levels(year),seq(2016,2026,5))))) %>%
  group_by(Iter,year) %>%
  summarise(x=sum(x)) %>%
  ungroup()

# Plotting #####################################################################
#

#Plots: internal migration#####
plot_od=ggplot(data=tod3111l, #%>% filter(Sex=="Males"),
       aes(x=year, y=x,  group=1,quantile=quantile)) +
  geom_fan(intervals=c(0.5,0.8,0.95)) +
  geom_interval(linetype="dotted", colour="darkgreen",intervals=c(0),show.legend = T) +
  geom_vline(xintercept = 2006, colour = "azure4", linetype="longdash") +
  geom_line(aes(x=year, y=dat_x, quantile=NULL, linetype="Data"), size=1.02,alpha=1,colour="black") + #colour="black"
  geom_line(data=filter(tod3111l,year>=2006),aes(x=year, y=dat_x, quantile=NULL), size=1.02,alpha=1,colour="blue",show.legend=F) + #, Sex=="Males"
  scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
  scale_linetype_manual(name="", values=c("solid","dotted")) +
  scale_y_continuous(name="Out-migration probability", #(common logarithm)
                     labels = function(x) format(round(10^(x),3), 
                                                 trim=T,digits=3)) +
  scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
  facet_grid2( Reg_O ~ Reg_D, scales = "free_y", switch = "y",independent = "y" ) + #,switch="y"
  # facet_grid( Reg_O ~ Reg_D, scales = "fixed", switch = "y") +
  guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0)))) +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold"), 
        axis.text.x  = element_text(size=11,angle=90, vjust=0.5), 
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.title = element_text(size=14), 
        strip.text.x = element_text(size=13), 
        strip.text.y = element_text(size=13),
        legend.position = "bottom")
ggsave("plots/Figure2.pdf",device="pdf",plot = plot_od,
       width = 6.4, height = 6.0, scale=1.6)#width = 4, height = 4, scale=2.5

#Plots: fertility####
plot_f=ggplot(data=tfr16,# %>% filter(S=="F"), #lex0.04 %>% select(-A)
       aes(x=year, y=x,  group=1,quantile=quantile)) +
  geom_fan(intervals=c(0.5,0.8,0.95)) +
  geom_interval(linetype="dotted", colour="darkgreen",intervals=c(0),size=.5,show.legend = T) +
  geom_line(aes(x=year, y=dat_x, quantile=NULL, linetype="Data"), size=1.02,alpha=1,colour="black",show.legend=F) + #, linetype="Data"
  geom_line(data=filter(tfr16,year>=2011),aes(x=year, y=dat_x, quantile=NULL), size=1.02,alpha=1,colour="blue",show.legend=F) +
  geom_vline(xintercept = 2011, colour = "azure4", linetype="longdash") +
  scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
  scale_size_manual(name="", values=c(1.02,0.5)) +
  scale_linetype_manual(name="", values=c("solid","dotted")) +
  scale_y_continuous(name="Total Fertility Rate") +
  scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
  facet_wrap( Reg ~., scales = "fixed", nrow=2) +
  guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
         linetype=guide_legend(override.aes = list(linetype = c("solid","dotted"), size=c(1.02,0.5)))) +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold"), 
        axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
        axis.text.y = element_text(size=13),
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.title = element_text(size=14), 
        strip.text.x = element_text(size=13), 
        strip.text.y = element_text(size=13),
        legend.position = "bottom"
# panel.grid.major.y = element_blank(),
)
ggsave("plots/Figure3.pdf",plot = plot_f, 
       device="pdf",width = 6.3, height = 5.3, scale=1.8)

#Plots: mortality####
plot.lex0<-function(inp=lex0.04,nam){
  ggplot(data=inp,# %>% filter(S=="F"), #lex0.04 %>% select(-A)
         aes(x=year, y=x, quantile=quantile)) +#
    geom_fan(aes(group=Sex),intervals=c(0.5,0.8,0.95)) +
    scale_fill_gradient(name="Estimates:\nPredictive Interval", low= c("#11BD11"), high=c("#BBFEBB"), breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
    geom_interval(linetype="dotted",colour="darkgreen",intervals=c(0),show.legend = T) +
    facet_wrap( Reg ~ ., scales = "free_y", nrow=2) +
    #scale_linetype_manual(values=c("twodash", "dotted"))
    # scale_linetype_manual() +
    geom_line(aes(x=year, y=dat_x, group=Sex, colour=Sex, quantile=NULL), size=1.02,alpha=1) +
    scale_linetype(name=c("")) +
    scale_colour_manual(name="Sex (data)", values=c("red","darkblue","blank")) +
    scale_y_continuous(name="Life expectancy at birth") +
    scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0)), order=2)) +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), 
          axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          legend.text = element_text(size=14),
          legend.title = element_text(size=14), 
          strip.text.x = element_text(size=13), 
          strip.text.y = element_text(size=13),
          legend.position = "bottom"
          ) 
}
plot_mort=plot.lex0(lex0.14,"02")
ggsave(paste0("Plots/Figure4_","02",".pdf"),plot = plot_mort,
       device="pdf",width = 6.3, height = 5.3, scale=1.8)
# characteristics for the median LEx
# filter(lex0.14,A=="0-4",year==2026,quantile%in%c(0.025,0.975), Sex=="Male",Reg=="ACT")


#Plots: immigration####
plot.imi.tot<-function(inp,nam){
  p=ggplot(data=inp,# %>% filter(S=="F"), #lex0.04 %>% select(-A)
           aes(x=year, y=x,  group=1,quantile=quantile)) +
    geom_fan(intervals=c(0.5,0.8,0.95)) +
    geom_interval(linetype="dotted", colour="darkgreen",intervals=c(0),show.legend = T) +
    geom_line(aes(x=year, y=dat_x, quantile=NULL, linetype="Data"), size=1.02,alpha=1,colour="black") + #colour="black"
    geom_line(data=filter(inp,year>=2011),aes(x=year, y=dat_x, quantile=NULL), size=1.02,alpha=1,colour="blue",show.legend=F) +
    geom_vline(xintercept = 2011, colour = "azure4", linetype="longdash") +
    scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
    scale_size_manual(name="", values=c(1.02,0.5)) +
    scale_linetype_manual(name="", values=c("solid","dotted")) +
    scale_y_continuous(name="immigration counts, 1000s",
                       labels = function(x)round(10^x/1000,0)) + #"Counts (common logarithm)"
    scale_x_continuous(name="", expand=expand_scale(mult = c(0, 0))) +
    facet_grid( Sex ~ Reg, scales = "fixed", switch = "y") +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
           linetype=guide_legend(override.aes = list(linetype = c("solid","dotted"), size=c(1.02,0.5)))) +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), 
          # axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=13),
          # axis.title.x = element_text(size=14), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14), 
          legend.text = element_text(size=14),
          legend.title = element_text(size=14), 
          strip.text.x = element_text(size=13), 
          strip.text.y = element_text(size=13),
          legend.position = "bottom") 
  return(p)
}
plot_imi=plot.imi.tot(timi031l,"031")
ggsave(paste0("plots/Figure5_","031",".pdf"),plot = plot_imi, device="pdf",width = 6.3, height = 4.5, scale=2)


#Plots: emigration####
plot_emi=ggplot(data=temi02l,# %>% filter(S=="F"), #lex0.04 %>% select(-A)
       aes(x=year, y=x,  group=1,quantile=quantile)) +
  geom_fan(intervals=c(0.5,0.8,0.95)) +
  geom_interval(linetype="dotted", colour="darkgreen",intervals=c(0),show.legend = T) +
  geom_line(aes(x=year, y=dat_x, quantile=NULL, linetype="Data"), size=1.02,alpha=1,colour="black") + #colour="black"
  geom_line(data=filter(temi02l,year>=2011),aes(x=year, y=dat_x, quantile=NULL), size=1.02,alpha=1,colour="blue",show.legend=F) +
  geom_vline(xintercept = 2011, colour = "azure4", linetype="longdash") +
  scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
  scale_size_manual(name="", values=c(1.02,0.5)) +
  scale_linetype_manual(name="", values=c("solid","dotted")) +
  scale_y_continuous(labels = function(x)round(10^x,2)) + #"Rate (common logarithm)"
  # annotation_logticks(scaled = T) +
  scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
  facet_grid( Sex ~ Reg, scales = "fixed", switch = "y") +
  guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
         linetype=guide_legend(override.aes = list(linetype = c("solid","dotted"), size=c(1.02,0.5)))) +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold"), 
        axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
        axis.text.y = element_text(size=13),
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.title = element_text(size=14), 
        strip.text.x = element_text(size=13), 
        strip.text.y = element_text(size=13),
        legend.position = "bottom") +
  labs(y="emigration rate")
# ggsave("plots/Figure6.pdf",device="pdf",width = 5, height = 3, scale=2.3)
g=ggarrange(plot_imi,NULL,plot_emi, ncol=1,nrow=3, labels = NULL, legend ="bottom", common.legend = TRUE, align = "hv",heights = c(1,-0.05,1))
ggexport(g,filename= paste0("Plots/Figure56_1.pdf"),width=14,height = 14)

# plots migration together####
ggplot(data=migrl10,# %>% filter(S=="F"), #lex0.04 %>% select(-A)
         aes(x=year, y=x,  group=1,quantile=quantile)) +
  geom_fan(intervals=c(0.5,0.8,0.95)) +
  geom_interval(linetype="dotted", colour="darkgreen",intervals=c(0),show.legend = T) +
  geom_line(aes(x=year, y=dat_x, quantile=NULL, linetype="Data"), size=1.02,alpha=1,colour="black") + #colour="black"
  geom_line(data=filter(migrl10,year>=2011),aes(x=year, y=dat_x, quantile=NULL), size=1.02,alpha=1,colour="blue",show.legend=F) +
  geom_vline(xintercept = 2011, colour = "azure4", linetype="longdash") +
  scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
  scale_size_manual(name="", values=c(1.02,0.5)) +
  scale_linetype_manual(name="", values=c("solid","dotted")) +
  scale_y_continuous(name="Common logarithm") + #"Rate (common logarithm)"
  scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
  facet_grid( component+Sex~ Reg, scales = "free_y",switch = "y") +
  guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
         linetype=guide_legend(override.aes = list(linetype = c("solid","dotted"), size=c(1.02,0.5)))) +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold"), 
        axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
        axis.text.y = element_text(size=13),
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.title = element_text(size=14), 
        strip.text.x = element_text(size=13), 
        strip.text.y = element_text(size=13),
        legend.position = "bottom") 
ggsave("plots/Figure56.pdf",device="pdf",width = 6.3, height = 3, scale=3)



#Plots: population by region/time####
plot.pop.time=function(nam){
  p = ggplot(data=Pf.df3,# %>% filter(S=="F"), #lex0.04 %>% select(-A)
             aes(x=year, y=x/100000,  group=1,quantile=quantile)) +
    geom_fan(intervals=c(0.5,0.8,0.95)) +
    geom_interval(linetype="dotted", colour="darkgreen",intervals=c(0),show.legend = T) +
    geom_line(data=ERP_pop30,aes(x=year, y=ERP/100000, quantile=NULL,linetype="Data"), size=1.02,alpha=1,colour="black") + #colour="black"
    scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
    geom_vline(xintercept = 2011, colour = "azure4", linetype="longdash") +
    scale_size_manual(name="", values=c(1.02,0.5)) +
    scale_linetype_manual(name="", values=c("solid","dotted")) +
    scale_y_continuous(name="Population in 100,000s") +
    scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
    facet_wrap2(Sex ~ Reg, nrow = 4,
                scales = "free_y", 
                strip = strip_split(c("right","top"))) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
           linetype=guide_legend(override.aes = list(linetype = c("solid","dotted"), size=c(1.02,0.5)))) +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), 
          axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          legend.text = element_text(size=14),
          legend.title = element_text(size=14), 
          strip.text.x = element_text(size=13), 
          strip.text.y = element_text(size=13),
          legend.position = "bottom",
          ) 
#  
}
# p_pop=plot.pop.time("v2.0.1")
plot_pop=plot.pop.time("v3.1")
# ggsave(paste0("plots/Figure7.pdf"),plot=plot_pop,device="pdf",width = 6.3, height = 2, scale=3)

#Plots: population by time total####
plot.pop.t.time=function(nam){
  p = ggplot(data=Pf.df4,# %>% filter(S=="F"), #lex0.04 %>% select(-A)
             aes(x=year, y=x/100000,  group=1,quantile=quantile)) +
    geom_fan(intervals=c(0.5,0.8,0.95)) +
    geom_interval(linetype="dotted", colour="darkgreen",intervals=c(0),show.legend = T) +
    geom_line(data=filter(ERP_pop4, year>1999),aes(x=year, y=ERP/100000, quantile=NULL,linetype="Data"), size=1.02,alpha=1,colour="black") + #colour="black"
    scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + #, minor_breaks=c(seq(0.5,0.95,0.1),0.95)
    geom_vline(xintercept = 2011, colour = "azure4", linetype="longdash") +
    scale_size_manual(name="", values=c(1.02,0.5)) +
    scale_linetype_manual(name="", values=c("solid","dotted")) +
    scale_y_continuous(name="") +
    scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
    # facet_wrap( Sex ~ Reg, scales = "free_y", nrow=2) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
           linetype=guide_legend(override.aes = list(linetype = c("solid","dotted"), size=c(1.02,0.5)))) +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), 
          axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          legend.text = element_text(size=14),
          legend.title = element_text(size=14), 
          strip.text.x = element_text(size=13), 
          strip.text.y = element_text(size=13),
          legend.position = "bottom") 
#   ggsave(paste0("plots/Figure8",nam,".pdf"),plot=p,device="pdf",width = 5, height = 3, scale=1.5)
}
plot_poptot=plot.pop.t.time("v3.1")
g21=ggarrange(plot_pop,plot_poptot, ncol=2,nrow=1, labels = "auto", legend ="none", common.legend = TRUE, widths = c(4,2.3))
ggexport(g21,filename= paste0("Plots/Figure7.pdf"),
         width=12,height = 9)


#Plots: Population by age/sex/region####
#plotting with quantiles
plot.pop.age=function(ye,nam){
  p = ggplot(data=Pf.df1 %>% filter(year%in%c(ye)),  #lex0.04 %>% select(-A) #%>% filter(year==ye)
             aes(x=Age, y=x/1000, group=Sex, color=Sex, quantile=quantile)) +
    geom_fan(aes(group=Sex), intervals=c(0.5,0.8,0.95)) +
    # geom_col(data=ERP_pop0 %>% filter(year==2016),aes(x=Age, y=ERP/1000, color=Sex), width=1, alpha=0) +
    # geom_line(data=ERP_pop0 %>% filter(year==2016),aes(x=Age, y=ERP/1000, linetype=Sex, quantile=NULL), size=1.02,alpha=1,colour="black") + #colour="black"
                #if (ye>2016) 2016 else ye)
    geom_line(aes(x=Age, y=ifelse(Sex=="Males",-ERP/1000,ERP/1000), linetype=Sex, quantile=NULL), size=1.02,alpha=1,colour="black") +
    scale_fill_gradient(name="Estimates\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) + 
    scale_linetype(name=paste0("ERP (",if (ye[1]>2016) 2016 else ye[1],")")) +
    scale_y_continuous(labels=function(x) abs(x), name="Population in 1000s") +
    scale_x_discrete(breaks=levels(Pf.df1$Age)[c(seq(1,18,3),18)],labels=levels(Pf.df1$Age)[c(seq(1,18,3),18)],expand=expand_scale(mult = c(0, 0))) +
    guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0)))) +
    # facet_grid( year ~ Reg, scales = "free_x") + 
    facet_wrap2(year ~ Reg, nrow = 4,
                scales = "free_x", 
                strip = strip_split(c("right","top"))) +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"), 
          axis.text.x  = element_text(size=13,angle=90, vjust=0.5), 
          axis.text.y = element_text(size=10),
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          legend.text = element_text(size=14),
          legend.title = element_text(size=14), 
          strip.text.x = element_text(size=13), 
          strip.text.y = element_text(size=13),
          legend.position = "bottom",
          panel.spacing.y = unit(1,"lines")) +
    coord_flip()
}  

plot_age1626=plot.pop.age(c(2016,2026),"v3.1")
ggsave(paste0("plots/Figure9_1626.pdf"),plot=plot_age1626,device="pdf",width = 5, height = 6, scale=2.2)

#unused
# g3=ggarrange(g21,p_pa1626, ncol=1,nrow=2, labels = c("","c"), legend ="bottom", common.legend = TRUE, heights = c(1.1,1))
# ggexport(g3,filename= paste0("Plots/Figure7_rta.pdf"),width=15.75,height = 9)


#Tables####
#internal migration forecasts description
#check sex
tod3111l %>% filter(Reg_O=="ACT",Reg_D=="NSW",quantile%in%c(0.025,0.5,0.975), year==2016) %>% View
tod3111l %>% filter(Reg_O=="NSW",Reg_D=="NT",quantile%in%c(0.025,0.5,0.975), year==2016) %>% View
#fertility TFR
tfr16 %>% filter(year==2016) %>% View
tfr16 %>% filter(year==2026) %>% View
#immigration
t_imi_tot=trans.imig.tot(t_imi44,FUN = identity)
t_imi_tot %>% filter(Reg=="NSW",year==2026)
#immigration totals plot
t_imitot %>% 
  ggplot(aes(x=year,y=x/1000,group=1,quantile=quantile)) + 
  geom_fan(intervals=c(0.5,0.8,0.95)) + 
  geom_interval(linetype="dotted", 
                colour="darkgreen",
                intervals=c(0),
                show.legend = T) +
  geom_vline(xintercept = 2011, colour = "azure4", linetype="longdash") +
  geom_line(aes(x=year,y=dat_x/1000,quantile=NULL,linetype="Data"), 
            size=1.02,alpha=1,colour="black") +
  scale_fill_gradient(name="Estimates:\nPredictive Interval", low= "#11BD11", high="#BBFEBB", breaks=c(0.5,0.8,0.95), labels=c("50%","80%","95%")) +
  scale_size_manual(name="", values=c(1.02,0.5)) +
  scale_linetype_manual(name="", values=c("solid","dotted")) +
  scale_y_continuous(name="immigration counts, 1000s") + 
  scale_x_continuous(name="", expand=expand_scale(mult = c(0, 0))) +
  guides(fill=guide_legend(override.aes = list(linetype = c(0,0,0))),
         linetype=guide_legend(override.aes = list(
           linetype = c("solid","dotted"), 
           size=c(1.02,0.5)))) +
  theme_bw() + 
  theme(plot.title = element_text(lineheight=.8, face="bold"), 
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size=13),
                        # axis.title.x = element_text(size=14), 
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size=14), 
                        legend.text = element_text(size=14),
                        legend.title = element_text(size=14), 
                        strip.text.x = element_text(size=13), 
                        strip.text.y = element_text(size=13),
                        legend.position = "bottom") 
ggsave(paste0("plots/Figure_imitot.pdf"),device="pdf",width = 6.3, height = 2.3, scale=2)

#ERP population totals by age/region/sex
filter(ERP_pop30, year==2011) %>% arrange(Sex)
filter(ERP_pop30, year==2016) %>% arrange(Sex)
#Table with population forecasts
Pf.df3 %>% spread(quantile, x) %>% filter(year==2016) %>% arrange(Sex) %>% select(-c(1:5,7,11)) %>% '/'(1000) %>% round(0)

PPD2016=Pf.df40 %>% filter(year==2016) %>% pull(x)
sum(PPD2016<=ERP_pop4 %>% filter(year==2016) %>% pull(ERP))/1000

