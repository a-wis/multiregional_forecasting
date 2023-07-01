#simulation

library(tidyverse)
library(rstan)

set.seed(26)
k1=c(0,1:19+rnorm(19,0,0.5))/20
k2=c(0,seq(0,-4,length.out=19)+rnorm(19,0,0.5))
a=c(-3,-2.5,-3,-2.3,-1.5)
b=c(-0.2,0.3,0.1,0.5,0.3)
r=c(0.5,0.1,.8,.9)-1
s=c(0.3,-0.3,0.6,0.4)

#ym=array(NA,dim=c(5,4,20))

AT=a+b%o%k1
RT=r+s%o%k2
yART=array(NA,dim = c(5,4,20),dimnames = list(A=c("5","10","15","20","25"),R=c("R1","R2","R3","R4"),Y=1:20))
for (i in 1:5){
  for (j in 1:4){
    for (k in 1:20){
      yART[i,j,k]=AT[i,k]+RT[j,k]+rnorm(1,0,0.1)
    }
  }
}
range(yART)
range(exp(yART))
dimnames(yART)
ydf=as.data.frame.array(yART) %>% tibble::rownames_to_column(var = "A") %>%
  pivot_longer(cols = 2:last_col(),names_to = c("R","Time"),names_sep = "\\.") %>%
  mutate(A=as_factor(A), Time=as.numeric(Time))

yA=ydf %>% group_by(A) %>%
  summarise(Am=mean(value)) %>%
  ungroup() %>% 
  arrange(A) 

yR=ydf %>%
  # left_join(yA) %>%
  # mutate(value=value-Am) %>%
  group_by(R) %>%
  summarise(Rm=mean(value)) %>%
  ungroup() %>% 
  arrange(R)

yT=ydf %>%
  group_by(Time) %>%
  summarise(Tm=mean(value)) %>%
  ungroup() %>% 
  arrange(Time)

prepA=ydf %>% 
  left_join(yR) %>%
  mutate(y=value-Rm) %>%
  group_by(A,Time) %>%
  summarise(value=mean(value)) %>%
  ungroup() %>% 
  left_join(yA) %>%
  mutate(y=value-Am) %>% 
  arrange(Time,A) %>% 
  select(-value,-Am) %>%
  pivot_wider(names_from = Time, values_from=y) %>%
  select(-A) %>%
  as.matrix()



foo1=prcomp(prepA, center = F, scale. = T)
foo1=prcomp(prepA, center = T, scale. = T)
plot(b,(foo1$x[,1])/sum((foo1$x[,1])))
# plot((foo1$x[,1])/sum((foo1$x[,1])))
# plot((foo1$x[,1]))
# plot(b,(foo1$x[,1])/sum((foo1$x[,1])))
# abline(a=0,b=1)

yAT=ydf %>% group_by(A,Time) %>%
  summarise(ATm=mean(value))
yRT=ydf %>% group_by(R,Time) %>%
  summarise(RTm=mean(value))

prepR=ydf %>% 
  left_join(yA) %>%
  mutate(y=value-Am) %>%
  group_by(R,Time) %>%
  summarise(value=mean(value)) %>%
  ungroup() %>% 
  left_join(yR) %>%
  mutate(y=value-Rm) %>% 
  arrange(Time,R) %>% 
  select(-value,-Rm) %>%
  pivot_wider(names_from = Time, values_from=y) %>%
  select(-R) %>%
  as.matrix()

# prepR=ydf %>% 
#   group_by(R,Time) %>%
#   summarise(value=mean(value)) %>%
#   ungroup() %>% 
#   left_join(yAT) %>%
#   mutate(value=value-ATm) %>%
#   arrange(Time,R) %>% 
#   pivot_wider(names_from = Time, values_from=value) %>%
#   select(-R) %>%
#   as.matrix()


foo2=prcomp(prepR,center = F,scale. =T)
foo2=prcomp(prepR,center = T,scale. =T)
plot(s,(foo2$x[,1])/sum((foo2$x[,1])))
# plot((foo2$x[,1])/sum((foo2$x[,1])))
# plot((foo2$x[,1]))
# plot(s,(foo2$x[,1])/sum((foo2$x[,1])))
# abline(a=0,b=1)


# modelling #####

data.inp.toy01<-list(y= ydf %>% pull(value),
                 mat_Rc=diag(3)+1, mat_Ra=diag(4)+1,
                 N=5,R=4,T=20,
                 sig_t=0.1)
#model 3_0 with OD_b
toyinits01 <- list(k1=c(1:19),k2=c(1:19), sig1=c(0.05),sig2=c(.5,.5,.5),sigk=c(0.1,0.1))
toyinits02 <- list(kn1=c(-10:9),kn2=c(-10:9), sig1=c(0.05),sig2=c(.5,.5),sigk=c(0.1,0.1))

system.time(fittoy_0161 <- stan(file = "models/stan_toy_0161.stan", data = data.inp.toy01, iter = 1500,verbose = FALSE, chains=2,thin=1,warmup = 1000,control = list(adapt_delta=0.98,max_treedepth=15),init = list(toyinits01,toyinits01),cores = 2,seed = 1349)
)

data.inp.toy02<-list(y= ydf %>% pull(value),
                     mat_Rc=diag(3)+1, mat_Ra=diag(4)+1,
                     A_b=(foo1$x[,1])/sum((foo1$x[,1])),
                     R_b=(foo2$x[,1])/sum((foo2$x[,1])),
                     # A_b=(foo1$x[,1]),
                     # R_b=(foo2$x[,1]),
                     # A_bp=(foo1$x[,1])/sum((foo1$x[,1])),
                     # R_bp=(foo2$x[,1])/sum((foo2$x[,1])),
                     k1_0=foo1$rotation[1,1],
                     k2_0=foo2$rotation[1,1],
                     N=5,R=4,T=20,
                     sig_t=5)
system.time(fittoy_0262 <- stan(file = "models/stan_toy_026.stan", data = data.inp.toy02, iter = 1500,verbose = FALSE, chains=2,thin=1,warmup = 1000,control = list(adapt_delta=0.95,max_treedepth=15),init = list(toyinits01,toyinits01),cores = 2,seed = 1349)
)

# COMMENTS ####
# model 016 performs relatively best, models with external PC
# have poorer fit ye they are faster, 
# adding a global constant seems to help?

plot(k2,c(0,summary(fittoy_016,pars="k2")$summary[,1]))
plot(k2,c(0,summary(fittoy_02,pars="k2")$summary[,1]))
plot(k1,c(0,summary(fittoy_016,pars="k1")$summary[,1]))
plot(k1,c(0,summary(fittoy_02,pars="k1")$summary[,1]))
plot(a,summary(fittoy_01,pars="A_a")$summary[,1])
plot(a,summary(fittoy_02,pars="A_a")$summary[,1])
plot(a,summary(fittoy_01,pars="R_a")$summary[,1])
plot(b,summary(fittoy_016,pars="A_b")$summary[,1])
plot(s,summary(fittoy_016,pars="R_b")$summary[,1])
plot(b,summary(fittoy_02,pars="A_b")$summary[,1])
plot(s,summary(fittoy_02,pars="R_b")$summary[,1])


plot(b,data.inp.toy02$A_b)
plot(s,data.inp.toy02$R_b)
plot(data.inp.toy02$A_b,summary(fittoy_01,pars="A_b")$summary[,1])
plot(data.inp.toy02$R_b,summary(fittoy_01,pars="R_b")$summary[,1])

plot(data.inp.toy01$y,summary(fittoy_016,pars="mmd")$summary[,1])
abline(0,1)

data.inp.toy02$A_b
summary(fittoy_02,pars="A_b")$summary[,1]
summary(fittoy_02,pars="k1")$summary[,1]

plot(fittoy_02,pars=c("k1"),plotfun="trace")
plot(fittoy_02,pars=c("k2"),plotfun="trace")
plot(fittoy_02,pars=c("k1"),plotfun="plot")
plot(fittoy_02,pars=c("k2"),plotfun="plot")
plot(fittoy_02,pars=c("k2"),plotfun="dens")
plot(fittoy_02,pars=c("sigk"),plotfun="trace")
plot(fittoy_02,pars=c("sigk"),plotfun="plot")
sum02=summary(fittoy_02)$summary
sum02[sum02[,10]>1.05,c(1,10)]

plot(fittoy_01,pars=c("A_a"),plotfun="dens")
plot(fittoy_02,pars=c("R_a"),plotfun="dens")
plot(fittoy_02,pars=c("A_a"),plotfun="dens")
plot(fittoy_02,pars=c("A_a"),plotfun="plot")
plot(fittoy_02,pars=c("R_a"),plotfun="plot")
plot(fittoy_02,pars=c("yf"),plotfun="plot")
plot(fittoy_02,pars=c("R_b"),plotfun="dens")
plot(fittoy_02,pars=c("A_b"),plotfun="dens")
plot(fittoy_02,pars=c("R_b"),plotfun="plot")
plot(fittoy_02,pars=c("A_b"),plotfun="plot")
