# esimation of log-bilinear models for population components
#A unified framework for probabilistic forecasting subnational populations
#(C) anonymised (2020) 

# preparation ####
library(rstan)
library(tidyr)
library(dplyr)
# library(shinystan)

#rounding function
round2=function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


#internal migration ####
#data
data.inp05<-list(d=input_au %>% filter(year!=2011) %>% pull(Flow),
                 p1=input_au %>% filter(year!=2011) %>% pull(x), 
                 mat_Rc=diag(55)+1, mat_Ra=diag(16)+1,
                 N=17,R=8,T=6,F=3,Co=56,S=2,
                 ind_DA=which(corridor!=8),
                 ind_DA0=which(corridor==8),
                 corridor=corridor,unicorr=unif_ind,
                 sig_t=0.5)
#inits
lineinits040 <- list(k1=c(1:5),k2=c(1:5), sig1=0.1, sig2=c(0.5,0.5),
                     sig3=c(.1,.1,.1),sigk=c(0.1,0.1))

#model - age/OD changing over time
system.time(fit08_04 <- stan(
  file = "models/stan_migr_internal_n8_04.stan", 
  data = data.inp05, iter = 1000,verbose = FALSE, 
  chains=2,thin=1, warmup = 500,
  control = list(adapt_delta=0.97,max_treedepth=18),
  init = list(lineinits040,lineinits040),cores = 2,seed = 1349)
)

#mortality ####
#data
data.inp_m11<-list(d=round2(deaths9311$x,0),p1=deaths9311$ERP, 
                   mat_RA=diag(17)+1, 
                   mat_RR=diag(7)+1,N=18,R=8,T=19,F=15,S=2, sig_t=0.5)

#inits
lineinits_mort_24 <- list(ben=matrix(c(rep(1:2,8),1)/sum(c(rep(1:2,9))),17,2), k=matrix(1:18,18,2), sig1=c(0.05),sig2=c(.5),sig3=c(.1,.1))

#model
system.time(fit_m02 <- stan(file = "models/stan_mort_02.stan", data = data.inp_m11, iter = 1500,verbose = FALSE, chains=2,thin=1,warmup = 1000,control = list(adapt_delta=0.96,max_treedepth=18),init = list(lineinits_mort_24,lineinits_mort_24),cores = 2,seed = 1349)
)


#emigration ####
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


#immigration ####
#data
data.inp_im01<-list(d=filter(intmigr_8116, 
                             year<2012)$immig, 
                    mat_RA=diag(17)+1, 
                    mat_RR=diag(7)+1, N=18,R=8,T=31,F=15,S=2, sig_t=0.5)

#inits
lineinits_imi_44 <- list(ben=matrix(c(rep(1:2,8),1)/sum(c(rep(1:2,9))),17,2), k=matrix(1:30,30,2), sig1=c(0.05),sig2=c(.5,.5,0.5),sig3=c(.1,.1))

#model
system.time(fit_imi031 <- stan(file = "models/stan_imig_031.stan", data = data.inp_im01, iter = 1300,verbose = FALSE, chains=2,thin=1,warmup = 800,control = list(adapt_delta=0.97,max_treedepth=15),init = list(lineinits_imi_4402, lineinits_imi_4402),cores = 2,seed = 1349)
)

#fertility ####
#data
data.inp_f01<-list(d=round(filter(births_8116,year<2012)$births),p1=filter(births_8116,year<2012)$ERP, mat_R=diag(6)+1,
                   mat_RR=diag(7)+1, 
                   N=7,R=8,T=31,F=15, sig_t=0.5)

#inits
lineinits_f_12 <- list(ben=c(rep(1:2,3))/sum(c(rep(1:2,3),1)), k=1:30,  sig1=c(0.05),sig2=c(0.5,0.5),sig3=c(0.1),coh=rep(0.1,50))

#model
system.time(fit_f01 <- stan(file = "models/stan_fert_01.stan", data = data.inp_f01, iter = 1800,verbose = FALSE, chains=2,thin=2,warmup = 1300,control = list(adapt_delta=0.99,max_treedepth=20),init = list(lineinits_f_12,lineinits_f_12),cores = 2,seed = 1349)
)

# convergence checks ####
plot(fit_f0016, plotfun = "rhat")
plot(fit_f0016, plotfun = "rhat", pars = c("mmd","mmd_f"))
plot(fit_f0016, plotfun = "trace", pars = c("sig1","sig2","sig3","sigk","sigc"), inc_warmup = T)
