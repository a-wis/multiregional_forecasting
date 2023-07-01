#libraries ####
library(broom)
library(HelpersMG)
#prep data for model selection ####
#year as factor

# RMSE function ####
rmse <- function(...,sort=FALSE) {
  calc_rmse <- function(model){
  fitted=model$fitted.values
  data_rate=model$data$x
  rmse=sqrt(mean((fitted-data_rate)^2))
  # mae=mean(abs(fitted-data_rate))
  # maxi=max(abs(fitted-data_rate))
  }
  
  # get_name<- function(...) as.character(as.list(match.call()[-1L]))
  
  result <- data.frame(model=as.character(as.list(match.call()[-1L])), 
             rmse=unlist(lapply(list(...), calc_rmse))) 
  if (sort==TRUE) result %>% dplyr::arrange(rmse)
  return(result)
}

input_au1 = input_au %>%
  mutate(year=as.factor(year))

## some test models####
mod1=glm(round2(Flow,0)~A+S+Corridor+year+A*Reg_O+A*Reg_D,family = poisson(link = "log"),data = input_au1,offset = log(x))

mod2=glm(round2(Flow,0)~A+S+Corridor+year+A*Reg_O+A*Reg_D+S*Reg_O+S*Reg_D,family = poisson(link = "log"),data = input_au,offset = log(x))

mod3=glm(round2(Flow,0)~A+S+Corridor+year+A*Reg_O+A*Reg_D+S*Reg_O+S*Reg_D + A*year,
         family = poisson(link = "log"),data = input_au1,
         offset = log(x))

##comparing LL models ####
### internal migration model selection ####
mod30=glm(round2(Flow,0)~(A+S+Corridor+year)^2,
         family = poisson(link = "log"),data = input_au1,
         offset = log(x))
mod31=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + Corridor*year + S*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod32=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + Corridor*year ),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod33=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + S*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod34=glm(round2(Flow,0)~(A+S+Corridor+year + A*year + Corridor*year + S*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod35=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod36=glm(round2(Flow,0)~(A+S+Corridor+year + A*year + Corridor*year ),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod37=glm(round2(Flow,0)~(A+S+Corridor+year + A*year + S*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod38=glm(round2(Flow,0)~(A+S+Corridor+year + A*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod39=glm(round2(Flow,0)~(A+S+Corridor+year + Corridor*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod40=glm(round2(Flow,0)~(A+S+Corridor+year + S*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))

#Corridor*A has perfect prediction (0s flat for some observations)
mod41=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + Corridor*A),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod42=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + Corridor*A + S*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod43=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + Corridor*A + Corridor*S),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod44=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + Corridor*A + Corridor*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
mod45=glm(round2(Flow,0)~(A+S+Corridor+year + A*S + A*year + Corridor*S + Corridor*year),
          family = poisson(link = "log"),data = input_au1,
          offset = log(x))
#modelling with O and D (not Corridor)
mod501=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D+year),
         family = poisson(link = "log"),data = input_au1,
         offset = log(x))
mod502=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D+ A*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod503=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D+ A*year + Corridor*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod504=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A ),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod505=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + Corridor*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod506=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*S + Reg_D*S + Corridor*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod507=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + Corridor*year + Reg_O*S + Reg_D*S),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod508=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + Corridor*year + Reg_O*S + Reg_D*S + S*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod509=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + Corridor*year + S*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod510=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + Corridor*year + A*S),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod511=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + Corridor*year + A*S + S*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod512=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + Corridor*year + A*S + Reg_O*S + Reg_D*S + S*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod513=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A +
                             Corridor*year + A*S + Reg_O*S + Reg_D*S),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod514=glm(round2(Flow,0)~(A+S+Reg_O+Reg_D + A*year + Reg_O*A + Reg_D*A + 
                             + A*S + Reg_O*S + Reg_D*S),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod515=glm(round2(Flow,0)~(A+S + A*year + Reg_O*A + Reg_D*A + Corridor*year + A*S*Reg_O + A*S*Reg_D ),#+ Reg_O*S + Reg_D*S - not needed
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod516=glm(round2(Flow,0)~(A+S + A*S*year + Reg_O*A + Reg_D*A + Corridor*year + Reg_O*S + Reg_D*S ),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
#equivalent to mod516
# mod517=glm(round2(Flow,0)~(A+S + A*S*year + Reg_O*A + Reg_D*A + Corridor*year + Reg_O*S + Reg_D*S + A*S),
           # family = poisson(link = "log"),data = input_au1,
           # offset = log(x))
# A*S*Corridor - zero counts - not identifiable

HelpersMG::compare_BIC(m1=mod501, m2=mod502, m3=mod503, 
                       m4=mod504, m5=mod505, m6=mod506, 
                       m7=mod507, m8=mod508, m9=mod509,
                       m10=mod510, m11=mod511, m12=mod512, 
                       m13=mod513, m14=mod514, m15=mod515,
                       m16=mod516)
rmse(mod501, mod502, mod503, 
     mod504, mod505, mod506, 
     mod507, mod508, mod509,
     mod510, mod511, mod512, 
     mod513, mod514, mod515,
     mod516)

#older
mod504=glm(round2(Flow,0)~(A+S+Corridor+year + A*year + Reg_O*A + Reg_D*A + Corridor*year+ A*S ),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod505=glm(round2(Flow,0)~(A+S+Corridor+year + A*year + Reg_O*A + Reg_D*A + Corridor*year+ S*year ),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod506=glm(round2(Flow,0)~(A+S+Corridor+year + A*year + Reg_O*A + Reg_D*A + Corridor*year+ S*Reg_O+ S*Reg_D ),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
#A+S+OD + OA + OD + AT + ODT
mod507=glm(round2(Flow,0)~(A+S+Corridor+year + A*year + Reg_O*A + Reg_D*A + A*S ),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))
mod5041=glm(round2(Flow,0)~(A+S+Corridor+ A*year + Reg_O*A + Reg_D*A + Corridor*year+ A*S ),
           family = poisson(link = "log"),data = input_au1,
           offset = log(x))

## information criteria ####
aic1=AIC(mod1, mod2, mod3, mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43, mod44,mod45, mod502,mod503,mod504,mod505,mod506)
bic1=BIC(mod1, mod2, mod3, mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43, mod44,mod45, mod502,mod503,mod504,mod505,mod506,mod507, mod5041)
View(bic1)
hqc1=AIC(mod1, mod2, mod3, mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43, mod44,mod45, mod502,mod503,mod504,mod505,mod506, k=2*log(log(nobs(mod504))))
# AICc = function(...){
#   AIC(..., k=2*(...$rank+1)/(nobs(...)-...$rank-1))
# }

mods=list(mod1, mod2, mod3, mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43, mod44,mod45, mod502,mod503,mod504,mod505,mod506)
aicc1=NULL
for (i in 1:24){
  aicc1=rbind(aicc1, AIC(mods[[i]],k=2*(mods[[i]]$rank+1)/(nobs(mods[[i]])-mods[[i]]$rank-1)) )
}
  View(cbind(aic1,bic1,hqc1,aicc1))

### mortality model selection ####

deaths9311 = deaths9311 %>%
    mutate(year=as_factor(year),
           mx=x/ERP)
deaths9311 %>% group_by(year) %>%
  summarise(mx=sum(x)/sum(ERP)) %>%
  ungroup() %>%  ggplot(aes(x=year,y=mx)) + geom_line(aes(group=1))

deaths9311 %>% group_by(A,S) %>%
  summarise(mx=log10(sum(x)/sum(ERP))) %>%
  ungroup() %>%  ggplot(aes(x=A,y=mx)) + geom_line(aes(colour=S,group=S))

deaths9311 %>% 
  mutate(state_name=as_factor(state_name)) %>%
  group_by(state_name,year,S) %>%
  summarise(mx=log10(sum(x)/sum(ERP))) %>%
  ungroup() %>% 
  group_by(state_name,S) %>%
  mutate(mxd=mx-lag(mx)) %>% 
  ggplot(aes(x=year,y=mx)) +
  # geom_line(aes(colour=Reg,group=Reg)) +
  # geom_line() +
  geom_line(aes(colour=S,group=S)) +
  facet_wrap(state_name~., nrow = 2)


# deaths9311 = deaths9311 %>%
  # mutate(year=as.numeric(year))

#loglin models
mod001=glm(round2(x,0)~(A+S+Reg+year),
          family = poisson(link = "log"),data = deaths9311,
          offset = log(ERP))
# tidy(mod001)   #model details
# glance(mod001) #coeffs
# augment(mod001) #data

# mod0011=glm(round2(x,0)~(A+S+Reg+as.numeric(year)),
#            family = poisson(link = "log"),data = deaths9311,
#            offset = log(ERP))
# coefplot::coefplot(mod001,vertical=T)

# mod002=glm(round2(x,0)~(A+S+Reg+year)^2,
#            family = poisson(link = "log"),data = deaths9311,
#            offset = log(ERP))
# coefplot::coefplot(mod002,vertical=T, intercept=F, pointSize=0)

mod002=glm(round2(x,0)~A+S+Reg+year + A*S + A*Reg + S*Reg + A*year + Reg*year + S*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
# summary(mod003)

mod003=glm(round2(x,0)~A+S+Reg+year + A*S + A*Reg + S*Reg + A*year + S*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))

mod004=glm(round2(x,0)~A+S+Reg+year + A*S + A*Reg + S*Reg + A*year + Reg*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))

mod005=glm(round2(x,0)~A+S+Reg+year + A*S + A*Reg + S*Reg + A*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))


mod006=glm(round2(x,0)~A+S+Reg+year + A*S + A*Reg +  A*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))

mod007=glm(round2(x,0)~A+S+Reg+year + A*S + A*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))

mod008=glm(round2(x,0)~A+S+Reg+year + A*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
#three-way AST
mod009=glm(round2(x,0)~A+S+Reg+year + A*Reg + S*Reg + A*S*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod010=glm(round2(x,0)~A+S+Reg+year + A*Reg + S*Reg + A*S*year + Reg*year ,#A*S ++ S*year - equiv
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))

mod011=glm(round2(x,0)~A+S+Reg+year + A*Reg + S*Reg + A*S*year + Reg*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))

mod013=glm(round2(x,0)~A+S+Reg+year + A*Reg + S*Reg + A*S*year + Reg*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod014=glm(round2(x,0)~A+S+Reg+year + A*Reg + S*Reg + A*S*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod015=glm(round2(x,0)~A+S+Reg+year + A*Reg + S*Reg + A*S*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
#coefplot::coefplot(mod014,vertical=T, intercept=F, pointSize=0.1)
#year as numeric
mod015=glm(round2(x,0)~A+S+Reg+ A*S + A*Reg + S*Reg + A*S*as.numeric(year),
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod016=glm(round2(x,0)~A+S+Reg+ A*S + A*Reg + S*Reg + A*S*as.numeric(year) + S*as.numeric(year),
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod017=glm(round2(x,0)~A+S+Reg+ A*S + A*Reg + S*Reg + A*S*as.numeric(year) + Reg*as.numeric(year),
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))

aic_mort=AIC(mod001, mod002, mod003, mod004, mod005, mod006, mod007, mod008, mod009, 
             mod010, mod011, mod012, mod013, mod014, mod015, mod016, mod017)
bic_mort=BIC(mod001, mod002, mod003, mod004, mod005, mod006, mod007, mod008, mod009, 
             mod010, mod011, mod012, mod013, mod014, mod015, mod016, mod017)

HelpersMG::compare_BIC(m1=mod001, m2=mod002,m3=mod003, m4=mod004, m5=mod005, m6=mod006, m7=mod007, m8=mod008, m9=mod009, m10=mod010)
#BIC(mod001, mod002, mod003, mod004, mod005, mod006, mod007, mod008, mod009, 
#    mod010, mod011, mod012, mod013, mod014)
rmse(mod001, mod002, mod003, mod004, mod005, mod006, mod007, mod008, mod009, mod010)

#fertility model selection ####
births_8116 = births_8116 %>%
  mutate(year=as_factor(year)) %>% 
  rename(x=births)

mod101=glm(round2(x,0)~(A+Reg+year),
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))
mod102=glm(round2(x,0)~A+Reg+year + A*Reg + A*year + Reg*year,
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))
mod103=glm(round2(x,0)~A+Reg+year + A*Reg + A*year,
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))
mod104=glm(round2(x,0)~A+Reg+year + A*year + Reg*year,
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))
mod105=glm(round2(x,0)~A+Reg+year + A*year,
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))

HelpersMG::compare_BIC(m1=mod101, m2=mod102, m3=mod103, m4=mod104, m5=mod105)
rmse(mod101, mod102, mod103, mod104, mod105)


#### ### ### ### ### ### ### ### ### ### ### ### ### ### # # # # # # ### # # ### # # ###
# comparisons with old results and data ####
# mod2=glm(round2(Flow,0)~-1+Corridor+A+S+year,family = poisson(link = "log"),data = input_au %>%
#            filter(Reg_O%in%c("1","2","3"),Reg_D%in%c("1","2","3")),offset = log(x))

dat_check=input_au %>% 
  mutate(mig_rate=round2(Flow,0)/x) 
dat_check=dat_check%>% 
  bind_cols(fitm1=mod1$fitted.values/dat_check$x,
            fitm2=mod2$fitted.values/dat_check$x) %>% 
  mutate(m1err=mig_rate-fitm1,
         m2err=mig_rate-fitm2)


dat_check %>% group_by(S) %>% summarise(mean(abs(m1err)),mean(abs(m2err)))


temp1=log(x)-rowMeans(log(x)) #tested with emifAU (A:18 x T:33)
foo=prcomp(temp1,center = F)
(-(foo$x[,1])/sum(-(foo$x[,1])))


#loading old results
load("C:/Users/msassaw5/Dropbox (The University of Manchester)/Research/Multiregional/Rfiles/fit03111_new.RData")
OD_b=summary(fit03111,pars="OD_b")$summary
A_b = summary(fit03111,pars="A_b")$summary


prepA=dat_check %>% group_by(A,year) %>%
  summarise(Flow=sum(Flow),x=sum(x)) %>%
  ungroup() %>%
  mutate(lrate=log(round2(Flow,0)))

prepA=prepA %>% group_by(A) %>% mutate(lrate_diff=lrate-lag(lrate)) 

ggplot(data=prepA, aes(x=A, y=(lrate))) + geom_line(aes(colour=year,group=year)) + geom_point(aes(colour=year,group=year))

ggplot(data=prepA, aes(x=year, y=(lrate))) + geom_line(aes(colour=A,group=A)) + geom_smooth(method="loess") + geom_point(aes(colour=A,group=A))
ggplot(data=prepA, aes(x=year, y=(lrate_diff))) + geom_line(aes(colour=A,group=A)) + geom_smooth(method="loess") + geom_point(aes(colour=A,group=A))

prepY=dat_check %>% group_by(year,Reg_O,Reg_D,Corridor) %>%
  summarise(Flow=sum(Flow),x=sum(x)) %>%
  ungroup() %>%
  mutate(lrate=log(round2(Flow,0))) %>%
  group_by(Corridor) %>%
  mutate(lratediff=lrate-lag(lrate)) %>%
  ungroup()
prepY %>% ggplot(aes(x=year,y=lratediff)) + geom_line() + facet_grid(Reg_O~Reg_D)
prepY %>% ggplot(aes(x=year,y=lrate)) + geom_line() + facet_grid(Reg_O~Reg_D)
prepY %>% ggplot(aes(x=year,y=lrate)) + geom_line(aes(colour=Corridor,group=Corridor)) + geom_smooth(aes(x=year,y=lrate),method="lm")#aes(colour=Corridor,group=Corridor)

prepA=prepA %>% 
  group_by(A) %>%
  summarise(lrateA=mean(lrate)) %>%
  right_join(prepA) %>%
  mutate(lrate0=lrate-lrateA) %>%
  arrange(year,A) %>%
  select(-Flow,-x,-lrateA,-lrate) %>%
  # select(-Flow,-x) %>%
  pivot_wider(names_from = year, values_from=lrate0) %>%
  select(-A) %>%
  as.matrix()

foo1=prcomp(prepA,center = F)
plot((-foo1$x[,1])/sum((foo1$x[,1])))
plot((foo1$x[,1]))
plot(A_b[,1],foo1$x[,1])
plot(fit03111,pars="A_b", plotfun="plot") + scale_y_reverse() + coord_flip()

CorrOD=dat_check %>% select(Reg_O,Reg_D,Corridor) %>% distinct()
prepOD=dat_check %>% group_by(Corridor,year) %>%
  summarise(Flow=sum(Flow),x=sum(x),mig_rate=mean(mig_rate)) %>%
  ungroup() %>%
  mutate(lrate=log(round2(Flow,0)/x)) %>%
  right_join(CorrOD)
# prepO=dat_check %>% group_by(Reg_O) %>%
#   summarise(Flow=sum(Flow),x=sum(x)) %>%
#   ungroup() %>%
#   mutate(lrateO=log(round2(Flow,0)/x)) %>%
#   right_join(CorrOD) %>%
#   select(-Flow,-x)
# prepD=dat_check %>% group_by(Reg_D) %>%
#   summarise(Flow=sum(Flow),x=sum(x)) %>%
#   ungroup() %>%
#   mutate(lrateD=log(round2(Flow,0)/x)) %>%
#   right_join(CorrOD) %>%
#   select(-Flow,-x)
prepOD=prepOD %>% 
  group_by(Corridor) %>%
  summarise(lrateOD=mean(lrate),
            lratemOD=mean(lratem)) %>%
  ungroup() %>%
  right_join(prepOD) %>%
  # right_join(prepO) %>%
  # right_join(prepD) %>%
  mutate(lrate0=lrate-lrateOD,
         lratem0=lratem-lratemOD) %>% #
  arrange(year,Corridor) %>% 
  select(-Flow,-x,-lrateOD,-lrate,-Reg_D,-Reg_O) %>%
  pivot_wider(names_from = year, values_from=lrate0) %>%
  select(-Corridor) %>%
  as.matrix()

#calculating principal component to compare with OD_b from old results
foo2=prcomp(prepOD,center = F)
plot((foo2$x[,1])/sum((foo2$x[,1])))
plot(-(foo2$x[,1]))
plot(OD_b[,1],foo2$x[,1])
plot(OD_b[,1],(foo2$x[,1])/sum((foo2$x[,1])))
plot(fit03111,pars="OD_b", plotfun="plot") + scale_y_reverse() + coord_flip()
plot(fit03111,pars="k2", plotfun="plot") + scale_y_reverse() + coord_flip()

# testing models

library(withr, lib.loc="C:\\Rlibrary")
library(rstan)

# testing a model with constraints forced by division ####
#data
data.inp02<-list(d=input_au %>% filter(year!=2011) %>% pull(Flow),
                 p1=input_au %>% filter(year!=2011) %>% pull(x), 
                 N=17,R=8,T=6,F=3,Co=56,S=2,
                 corridor=corridor,unicorr=unif_ind,
                 ind_DA=which(corridor!=8),
                 ind_DA0=which(corridor==8))
lineinits032 <- list(k1=c(1:5),k2=c(1:5), sig1=c(0.05),sig2=c(.5,.5),sig3=c(.1,.1,.1,0.1),sigk=c(0.1,0.1))

system.time(fit033 <- stan(file = "models/stan_migr_internal_n3.stan", data = data.inp02, iter = 1700,verbose = FALSE, 
init = list(lineinits032,lineinits032),cores = 2,chains=2,thin=1,warmup = 1200,control = list(adapt_delta=0.98,max_treedepth=15),seed = 1349)
)

# testing a small model ####
input_test = input_au %>% filter(year!=2011, Reg_O%in%c("1","2","3"),Reg_D%in%c("1","2","3"), A%in%c("5","10","15","20","25","30"))
corridor=as.numeric(input_test$Reg_D[1:6])
#creating index for unidirectional flows
uniflow=as.numeric(input_test$uniflow)
unif_ind=1
for (i in 2:6) unif_ind[i]=if(prod(uniflow[i]!= uniflow[1:(i-1)])>0) max(unif_ind)+1 else unif_ind[which(uniflow[1:(i-1)]==uniflow[i])]
rm(i)

data.inp02<-list(d= input_test %>% pull(Flow),
                 p1=input_test %>% pull(x), 
                 N=6,R=3,T=6,F=3,Co=6,S=2,
                 corridor=corridor,unicorr=unif_ind,
                 ind_DA=which(corridor!=3),
                 ind_DA0=which(corridor==3))
#model 3_0 with OD_b
lineinits032 <- list(k1=c(1:5),k2=c(1:5), sig1=c(0.05),sig2=c(.5,.5),sig3=c(.1,.1,.1,0.1),sigk=c(0.1,0.1))
#model 3_1 without OD_b
lineinits032 <- list(k1=c(1:5),k2=c(1:5), sig1=c(0.05),sig2=c(.5,.5),sig3=c(.1,.1,.1,0.1),sigk=c(0.1))

system.time(fit033_0 <- stan(
  file = "models/stan_migr_internal_n3.stan", 
  data = data.inp02, iter = 1700,verbose = FALSE, 
  init = list(lineinits032,lineinits032),
  cores = 2,chains=2,thin=1,warmup = 1200,
  control = list(adapt_delta=0.98,max_treedepth=15),seed = 1349)
)

data.inp04<-list(d=input_test %>% pull(Flow),
                 p1=input_test %>% pull(x), 
                 mat_Rc=diag(5)+1, mat_Ra=diag(5)+1,
                 N=6,R=3,T=6,F=3,Co=6,S=2,
                 corridor=corridor,unicorr=unif_ind,
                 ind_DA=which(corridor!=3),
                 ind_DA0=which(corridor==3),
                 sig_t=5)
 

lineinits040 <- list(k1=c(1:5),k2=c(1:5), sig1=0.1, sig2=c(0.5),sig3=c(.1,.1,.1,0.1),sigk=c(0.1,0.1))

system.time(fit04_00 <- stan(
  file = "models/stan_migr_internal_n4_00.stan", 
  data = data.inp04, iter = 1700,verbose = FALSE, 
  chains=2,thin=1,warmup = 1200,
  control = list(adapt_delta=0.98,max_treedepth=15),
  init = list(lineinits040,lineinits040),cores = 2,seed = 1349)
)







# model with constraints with multivariate normal (Wisniowski et al. 2015)
data.inp02<-list(d=input_au %>% filter(year!=2011) %>% pull(Flow),
                 p1=input_au %>% filter(year!=2011) %>% pull(x), 
                 mat_R=diag(55)+1,N=17,R=8,T=6,F=3,Co=56,S=2,
                 corridor=corridor,unicorr=unif_ind)
flows_sum_age = input_au %>%
  group_by(A) %>%
  summarise_at(vars(Flow,x),sum) %>%
  mutate(rate=Flow/x) 

lineinits031 <- list(k1=c(1:5),k2=c(1:5), A_a=log(flows_sum_age$rate), sig1=c(0.05),sig2=c(.5,.5,0.5,0.5),sig3=c(.1,.1,.1,0.1),sigk=c(0.1,0.1))

system.time(fit03111n <- stan(
  file = "models/stan_migr_internal_n2.stan", 
  data = data.inp02, iter = 1700,verbose = FALSE, 
  chains=2,thin=1,warmup = 1200,
  control = list(adapt_delta=0.98,max_treedepth=15),
  init = list(lineinits031,lineinits031),cores = 2,seed = 1349)
)

# models with parsimoneous specs
data.inp04<-list(d=input_au %>% filter(year!=2011) %>% pull(Flow),
                 p1=input_au %>% filter(year!=2011) %>% pull(x), 
                 mat_Rc=diag(55)+1, mat_Ra=diag(16)+1,
                 N=17,R=8,T=6,F=3,Co=56,S=2,
                 ind_DA=which(corridor!=8),
                 ind_DA0=which(corridor==8),
                 corridor=corridor,unicorr=unif_ind,
                 sig_t=5)

lineinits040 <- list(k1=c(1:5),k2=c(1:5), sig1=0.1, sig2=c(0.5),sig3=c(.1,.1,.1,0.1),sigk=c(0.1,0.1))

system.time(fit04_01 <- stan(
  file = "models/stan_migr_internal_n4_01.stan", 
  data = data.inp04, iter = 1800,verbose = FALSE, 
  chains=2,thin=1, warmup = 1300,
  control = list(adapt_delta=0.95,max_treedepth=18),
  init = list(lineinits040,lineinits040),cores = 2,seed = 1349)
)

library(loo)
log_lik <- extract_log_lik(fit04_00, merge_chains = F)
r_eff <- relative_eff(exp(log_lik), cores = 4)
loo_0 <- loo(log_lik, r_eff =r_eff,cores = 4)
print(loo_0)

#testing t-distribution ####
data.t=list(nu=c(1, 2, 2.5, 3, 5), sigma=c(2,5))

t_dist=stan(
  file = "models/stan_t-dist.stan", 
  data = data.t, iter = 1000,verbose = FALSE, 
  chains=2,thin=1, warmup = 500,
  control = list(adapt_delta=0.95,max_treedepth=18),
  cores = 2,seed = 1349, algorithm = "Fixed_param")

#forecasting internal migration ####

input_au %>%  group_by(year,Corridor) %>% 
  summarise(asmr=log(sum(Flow)/sum(x))) %>%  
  group_by(Corridor) %>% 
  mutate(asmr_d=asmr-lag(asmr)) %>% 
  summarise(asmr_dm=mean(asmr_d, na.rm=T)) %>% 
  ggplot() + 
  geom_point(aes(x=Corridor,y=-asmr_dm,group=1), size=2)

#library(withr, lib.loc="C:\\Rlibrary")
#.libPaths(c("C:/Program Files/R/R-4.0.5/library",.libPaths()))
.libPaths(c("C:/Rlibrary/4.2.1",.libPaths()))
library(rstan)
options(mc.cores = 2)
rstan_options(disable_march_warning = TRUE)
# +SA interaction
data.inp05<-list(d=input_au %>% filter(year!=2011) %>% pull(Flow),
                 p1=input_au %>% filter(year!=2011) %>% pull(x), 
                 mat_Rc=diag(55)+1, mat_Ra=diag(16)+1,
                 N=17,R=8,T=6,F=3,Co=56,S=2,
                 ind_DA=which(corridor!=8),
                 ind_DA0=which(corridor==8),
                 corridor=corridor,unicorr=unif_ind,
                 sig_t=0.5)

lineinits040 <- list(k1=c(1:5),k2=c(1:5), sig1=0.1, sig2=c(0.5,0.5),sig3=c(.1,.1,.1),sigk=c(0.1,0.1))

system.time(fit08_02 <- stan(
  file = "models/stan_migr_internal_n8_02.stan", 
  data = data.inp05, iter = 1000,verbose = FALSE, 
  chains=2,thin=1, warmup = 500,
  control = list(adapt_delta=0.97,max_treedepth=18),
  init = list(lineinits040,lineinits040),cores = 2,seed = 1349)
)

log_lik <- extract_log_lik(fit05_01, merge_chains = F)
r_eff <- relative_eff(exp(log_lik), cores = 4)
loo_0 <- loo(log_lik, r_eff =r_eff,cores = 4)
print(loo_0)

#forecasting immigration ####
#data
data.inp_im01<-list(d=filter(intmigr_8116, 
                             year<2012)$immig, 
                    mat_RA=diag(17)+1, 
                    mat_RR=diag(7)+1, N=18,R=8,T=31,F=15,S=2, sig_t=5)

#inits
# lineinits_imi_44 <- list(ben1=matrix(c(rep(1:2,8),1)/sum(c(rep(1:2,9))),17,2),
#                          ben2=rep(1:7)/sum(1:7),
#                          k1=matrix(1:30,30,2), sig1=c(0.05),sig2=c(.5,.5,0.5,.5),sig3=c(.1,.1))

lineinits_imi_4402 <- list(ben=matrix(c(rep(1:2,8),1)/sum(c(rep(1:2,9))),17,2), k=matrix(1:30/10,30,2), sig1=c(0.2),sig2=c(5,.5,1,.1),sig3=c(.1,.1))


#model
system.time(fit_imi001 <- stan(file = "models/stan_imig_01.stan", data = data.inp_im01, iter = 1300,verbose = FALSE, chains=2,thin=1,warmup = 800,control = list(adapt_delta=0.97,max_treedepth=15),init = list(lineinits_imi_4402, lineinits_imi_4402),cores = 2,seed = 1349)
)


foo1=summary(fit_imi001,pars=c("k2"),plotfun="trace")
plot(foo1$c_summary[,2,1])
plot(foo1$c_summary[,2,2])
foo2=summary(fit_imi001,pars=c("R_b"),plotfun="trace")

dat=intmigr_8116 %>% filter(year<2012) %>% group_by(year,Reg) %>% summarise(immig=sum(immig)) %>% pivot_wider(names_from = Reg, values_from = immig, names_prefix = "R")
datm=as.matrix(log(dat))
datm_mean=colMeans(datm)
inp=(t(log(datm))-datm_mean)[-1,]
prc1=prcomp(inp,center=F)
plot(inp[1,],t="l",ylim=range(inp))
for (i in 2:8) lines(inp[i,],t="l",col=i)
plot(inp[,1],t="l",ylim=range(inp))
for (i in 2:31) lines(inp[,i],t="l",col=i)


#forecasting mortality ####
deaths9311 %>%  group_by(year,A) %>% summarise(asmr=log(sum(x)/sum(ERP))) %>%  group_by(A) %>% mutate(asmr_d=asmr-lag(asmr)) %>% summarise(asmr_dm=mean(asmr_d, na.rm=T)) %>% ggplot() + geom_line(aes(x=A,y=-asmr_dm,group=1), size=1.25)
deaths9311 %>%  group_by(year,Reg) %>% summarise(asmr=log(sum(x)/sum(ERP))) %>%  group_by(Reg) %>% mutate(asmr_d=asmr-lag(asmr)) %>% summarise(asmr_dm=mean(asmr_d, na.rm=T)) %>% ggplot() + geom_point(aes(x=Reg,y=-asmr_dm,group=1), size=2)
#data
data.inp_m11<-list(d=round2(deaths9311$x,0),p1=deaths9311$ERP, 
                   mat_RA=diag(17)+1, 
                   mat_RR=diag(7)+1,N=18,R=8,T=19,F=15,S=2, sig_t=0.5)

#inits
lineinits_mort_24 <- list(ben=matrix(c(rep(1:2,8),1)/sum(c(rep(1:2,9))),17,2), k=matrix(1:18,18,2), sig1=c(0.05),sig2=c(.5),sig3=c(.1,.1))

#model
system.time(fit_m02 <- stan(file = "models/stan_mort_02.stan", data = data.inp_m11, iter = 1500,verbose = FALSE, chains=2,thin=1,warmup = 1000,control = list(adapt_delta=0.96,max_treedepth=18),init = list(lineinits_mort_24,lineinits_mort_24),cores = 2,seed = 1349)
)
#forecasting emigration ####
#data
data.inp_em01<-list(d=filter(intmigr_8116, year<2012)$emig,
                    p1=filter(intmigr_8116,year<2012)$ERP, 
                    mat_RA=diag(17)+1, 
                    mat_RR=diag(7)+1,N=18,R=8,T=31,F=15,S=2,
                    sig_t=0.5)

#inits
lineinits_emi_44 <- list(ben=c(rep(1:2,8),1)/sum(c(rep(1:2,9))), k=matrix(1:30,30,2), sig1=c(0.2),sig2=c(.5,.5,.5),sig3=c(.1,.1))

#model
system.time(fit_emi03 <- stan(file = "Models/stan_emig_03.stan", data = data.inp_em01, iter = 1300,verbose = FALSE, chains=2,thin=1,warmup = 800,control = list(adapt_delta=0.95,max_treedepth=18),init = list(lineinits_emi_44,lineinits_emi_44),cores = 2,seed = 1349)
)

# forecasting fertility ####

births_8116 %>% group_by(year,A) %>% summarise(asfr=sum(births)/sum(ERP)) %>% ggplot() + geom_line(aes(x=year,y=asfr,colour=A), size=1.25) 
births_8116 %>% group_by(year,Reg) %>% summarise(asfr=sum(births)/sum(ERP)) %>% ggplot() + geom_line(aes(x=year,y=asfr,colour=Reg), size=1.25) 

births_8116 %>% group_by(year,A) %>% summarise(asfr=log(sum(births)/sum(ERP))) %>%  group_by(A) %>% mutate(asfr_d=asfr-lag(asfr)) %>% summarise(asfr_dm=mean(asfr_d,na.rm=T)) %>% ggplot() + geom_line(aes(x=A,y=asfr_dm,group=1), size=1.25) 
births_8116 %>% group_by(year,Reg) %>% summarise(asfr=log(sum(births)/sum(ERP))) %>%  group_by(Reg) %>% mutate(asfr_d=asfr-lag(asfr)) %>% summarise(asfr_dm=mean(asfr_d,na.rm=T)) %>% ggplot() + geom_point(aes(x=Reg,y=asfr_dm,group=1), size=2.25) 


intmigr_8116 %>% group_by(year,Reg) %>% summarise(asir=log(sum(immig))) %>%  group_by(Reg) %>% mutate(asir_d=asir-lag(asir)) %>% summarise(asir_dm=mean(asir_d,na.rm=T)) %>% ggplot() + geom_point(aes(x=Reg,y=asir_dm,group=1), size=2.25)
intmigr_8116 %>% group_by(year,Reg) %>% summarise(asir=log(sum(emig)/sum(ERP))) %>%  group_by(Reg) %>% mutate(asir_d=asir-lag(asir)) %>% summarise(asir_dm=mean(asir_d,na.rm=T)) %>% ggplot() + geom_point(aes(x=Reg,y=asir_dm,group=1), size=2.25)

#data
data.inp_f01<-list(d=round(filter(births_8116,year<2012)$births),p1=filter(births_8116,year<2012)$ERP, mat_R=diag(6)+1,mat_RR=diag(7)+1, N=7,R=8,T=31,F=15, sig_t=0.5)

#inits
lineinits_f_12 <- list(ben=c(rep(1:2,3))/sum(c(rep(1:2,3),1)), k=1:30,  sig1=c(0.05),sig2=c(0.5,0.5),sig3=c(0.1),coh=rep(0.1,50))

#model
system.time(fit_f01 <- stan(file = "models/stan_fert_01.stan", data = data.inp_f01, iter = 1800,verbose = FALSE, chains=2,thin=2,warmup = 800,control = list(adapt_delta=0.99,max_treedepth=20),init = list(lineinits_f_12,lineinits_f_12),cores = 2,seed = 1349)
)
