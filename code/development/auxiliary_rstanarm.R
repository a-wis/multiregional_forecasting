######### Log-linear rstanarm ####
# modelling using rstanarm to explore log-linear models without bilinear component
library(rstanarm)
library(bayesplot)
#year as factor
input_au= input_au %>% mutate(year= as.factor(year))

# main effects
mod1= stan_glm(Flow ~ Reg_O + Reg_D + A + S + year - 1,
               data = input_au %>% filter(year!=2011),
               offset = log(x),
               prior = normal(0,5),
               chains=2, cores=2,
               warmup=500, iter=1500,
               family = poisson(link="log"))
summary(mod1)


mcmc_areas(mod1)
mcmc_intervals(mod1)

# A interactions
mod2= stan_glm(Flow ~ Reg_O*A + Reg_D*A + A * S + A* year - 1,
               data = input_au %>% filter(year!=2011),
               offset = log(x),
               prior = normal(0,5),
               chains=2, cores=2,
               warmup=2000, iter=2500,
               family = poisson(link="log"))
save(mod2, file="outputs/mod2.RData")
# A + OD
mod21= stan_glm(Flow ~ Corridor + Reg_O*A + Reg_D*A + A * S + A* year - 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2300, iter=2800,
                family = poisson(link="log"))
save(mod21, file="outputs/mod21.RData")
c21=mcmc_intervals(mod21)
ggsave(filename = "outputs/mod21coef.pdf",plot = c21,device = "pdf", height = 35,width=6)
# S
mod22= stan_glm(Flow ~ Reg_O*S + Reg_D*S + A * S + A* year + S* year- 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2000, iter=2500,
                family = poisson(link="log"))
save(mod22, file="outputs/mod22.RData")

# S + OD
mod23= stan_glm(Flow ~ Reg_O*Reg_D + Reg_O*S + Reg_D*S + A*S + A*year + S*year - 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2000, iter=2500,
                family = poisson(link="log"))

# A+S -ST
mod24= stan_glm(Flow ~ Reg_O*A + Reg_D*A + Reg_O*S + Reg_D*S + A*S + A*year - 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2000, iter=2500,
                family = poisson(link="log"))
save(mod24, file="outputs/mod24.RData")
# A+S +ST
mod241= stan_glm(Flow ~ Reg_O*A + Reg_D*A + Reg_O*S + Reg_D*S + A*S + A*year + S*year- 1,
                 data = input_au %>% filter(year!=2011),
                 offset = log(x),
                 prior = normal(0,5),
                 chains=2, cores=2,
                 warmup=2000, iter=2500,
                 family = poisson(link="log"))
save(mod241, file="outputs/mod241.RData")

mod242= stan_glm(Flow ~ Corridor + Reg_O*A + Reg_D*A + Reg_O*S + Reg_D*S + A*S + A*year - 1,
                 data = input_au %>% filter(year!=2011),
                 offset = log(x),
                 prior = normal(0,5),
                 chains=2, cores=2,
                 warmup=2500, iter=3000,
                 
                 family = poisson(link="log"))
save(mod24, file="outputs/mod24.RData")


pp_check(mod1,plotfun = "stat_2d") + xlim(c(385,390))
pp_check(mod2,plotfun = "stat_2d") + xlim(c(385,390))
pp_check(mod21,plotfun = "stat_2d") + xlim(c(385,390)) + ylim(c(690,750))
pp_check(mod22,plotfun = "stat_2d") + xlim(c(385,390))+ ylim(c(690,750))
pp_check(mod24,plotfun = "stat_2d") + xlim(c(385,390))
pp_check(mod241,plotfun = "stat_2d") + xlim(c(385,390))+ ylim(c(690,750))
pp_check(mod32,plotfun = "stat_2d") + xlim(c(385,390)) + ylim(c(690,750))
library(loo)
loo1=loo(mod1)
loo2=loo(mod2)
loo21=loo(mod21)
loo22=loo(mod22)
loo23=loo(mod23)
loo24=loo(mod24)
loo241=loo(mod241)
loo32=loo(mod32)

# A+S + OD - ST
mod25= stan_glm(Flow ~ Corridor + Reg_O*A + Reg_D*A + Reg_O*S + Reg_D*S + A*S + A* year - 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=3000, iter=4000,
                control = list(adapt_delta=0.9,max_treedepth=20),
                family = poisson(link="log"))
save(mod25, file="outputs/mod25.RData")
# A+S + OD + ST
mod251= stan_glm(Flow ~ Corridor + Reg_O*A + Reg_D*A + Reg_O*S + Reg_D*S + A*S + A*year + S*year - 1,
                 data = input_au %>% filter(year!=2011),
                 offset = log(x),
                 prior = normal(0,5),
                 chains=2, cores=2,
                 warmup=2500, iter=3000,
                 control = list(adapt_delta=0.9,max_treedepth=15),
                 family = poisson(link="log"))
save(mod251, file="outputs/mod251.RData")

# A+S + OD + ST
mod252= stan_glm(Flow ~ Reg_O*Reg_D + Reg_O*A + Reg_D*A + Reg_O*S + Reg_D*S + A*S + A* year - 1,
                 data = input_au %>% filter(year!=2011),
                 offset = log(x),
                 prior = normal(0,5),
                 chains=2, cores=2,
                 warmup=2000, iter=2500,
                 family = poisson(link="log")) 

summary(mod2)

mcmc_intervals(mod2)
mcmc_areas(mod2,regex_pars = "year")
mcmc_intervals(mod2,regex_pars = "year")

mod21= stan_glm(Flow ~ Reg_O + Reg_D + A * S + A* year + S* year- 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=1000, iter=1500,
                family = poisson(link="log"))

mod22= stan_glm(Flow ~ Reg_O*year + Reg_D*year + A * S + A* year + S* year- 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2000, iter=2500,
                family = poisson(link="log"))

mod22=glm(Flow ~ Reg_O * Reg_D + A * S + A* year + S* year- 1,
          data = input_au %>% filter(year!=2011),
          offset = log(x),
          # prior = normal(0,5),
          # chains=2, cores=2,
          # warmup=1000, iter=1500,
          family = poisson(link="log"))


mod22=glm(Flow ~ Reg_O * Reg_D + A * S + A* year + S* year- 1,
          data = input_au %>% filter(year!=2011),
          offset = log(x),
          # prior = normal(0,5),
          # chains=2, cores=2,
          # warmup=1000, iter=1500,
          family = poisson(link="log"))

mod31= stan_glm(Flow ~ Reg_O*year + Reg_D*year + A * S * year - 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2000, iter=2500,
                family = poisson(link="log"))

mod31= glm(Flow ~ Reg_O*year + Reg_D*year + A * S * year - 1,
           data = input_au %>% filter(year!=2011),
           offset = log(x),
           # prior = normal(0,5),
           # chains=2, cores=2,
           # warmup=2000, iter=2500,
           family = poisson(link="log"))

mod32= stan_glm(Flow ~ Reg_O*Reg_D + A * S * year - 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2000, iter=2500,
                family = poisson(link="log"))

mod32= stan_glm(Flow ~ Corridor + A * S * year - 1,
                data = input_au %>% filter(year!=2011),
                offset = log(x),
                prior = normal(0,5),
                chains=2, cores=2,
                warmup=2000, iter=2500,
                family = poisson(link="log"))


pp_check(mod31,plotfun = "stat_2d") + xlim(c(385,390))

mod4= glm(Flow ~ Reg_O*Reg_D*year*A * S - 1,
          data = input_au %>% filter(year!=2011),
          offset = log(x),
          # prior = normal(0,5),
          # chains=2, cores=2,
          # warmup=2000, iter=2500,
          family = poisson(link="log"))



