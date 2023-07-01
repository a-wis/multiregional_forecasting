#libraries ####
library(broom)
library(HelpersMG)

#prep data for model selection ####
#year as factor

# functions ####
## RMSE function ####
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

##CV function ####
cv_rmse <- function(..., prop=0.9,n=5,data){
  mods=list(...)
  errs=matrix(NA,length(mods),n)
  for (j in 1:n){
    train=data %>% slice_sample(prop=prop)
    test=anti_join(data,train)
    rmse=NULL
    for (i in 1:length(mods)){
      mod=glm(formula = mods[[i]]$formula,
              family = poisson(link = "log"),data = train,
              offset = log(ERP))
      preds=predict.glm(object = mod, newdata = test,type = "response")
      rmse=c(rmse,sqrt(sum((preds-test$x)^2)))
    }
    errs[,j] <- rmse
  }
  result <- data.frame(model=head(as.character(as.list(match.call()[-1L])),-3),
                       cv_rmse=rowMeans(errs))
  #
  return(result)
}

input_au1 = input_au %>%
  mutate(year=as.factor(year)) %>%
  rename(ERP=x, x=Flow)

deaths9311 = deaths9311 %>%
  mutate(year=as_factor(year),
         mx=x/ERP)

#internal migration model selection ####
mod501=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
mod502=glm(round2(x,0)~(A+S+Reg_O+Reg_D+ A*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
mod503=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
#CHOOSE
mod504=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Reg_O*A+ Reg_D*A),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
mod505=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Reg_O*S+ Reg_D*S),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
mod506=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Reg_O*S+ Reg_D*S+ Reg_O*A+ Reg_D*A),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
#AST
mod507=glm(round2(x,0)~(A+S+Reg_O+Reg_D+ Reg_O*S+ Reg_D*S+ Reg_O*A+ Reg_D*A+A*S*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
#504+ODT
mod508=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Reg_O*A+ Reg_D*A + Corridor*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
#504+ST
mod509=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Reg_O*A+ Reg_D*A + S*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
#504 +RT+ST
mod510=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Reg_O*A+ Reg_D*A + S*year+ Corridor*year),
           family = poisson(link = "log"),data = input_au1,
           offset = log(ERP))
#ODA
# mod511=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Corridor*A +  Corridor*year),
           # family = poisson(link = "log"),data = input_au1,
           # offset = log(ERP))
# mod512=glm(round2(x,0)~(A+S+Reg_O+Reg_D+year+ A*year + A*S + Corridor*S +  Corridor*year),
#            family = poisson(link = "log"),data = input_au1,
#            offset = log(ERP))

rmse(mod501,mod502,mod503,mod504,mod505,mod506,mod507,mod508,mod509,mod510)#,mod511)
HelpersMG::compare_BIC(m1=mod501,m2=mod502,m3=mod503,m4=mod504,m5=mod505,m6=mod506,m7=mod507,m8=mod508,m9=mod509,m10=mod510)

# mortality model selection ####
#loglin models
deaths9311 = deaths9311 %>%
  mutate(year=as_factor(year)) 
  
mod001=glm(round2(x,0)~(A+S+Reg+year),
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
# tidy(mod001)   #model details
# glance(mod001) #coeffs
# augment(mod001) #data
#AT
mod002=glm(round2(x,0)~A+S+Reg+year + A*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
#AS 
mod003=glm(round2(x,0)~A+S+Reg+year + A*year + A*S,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
#R effects
mod004=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod005=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + S*Reg,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod006=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg + S*Reg,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
#rmse() mod004
mod007=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod008=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + S*Reg,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))
mod009=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + S*Reg + Reg*year,
           family = poisson(link = "log"),data = deaths9311,
           offset = log(ERP))


#S*year is already in m008
# mod008=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + Reg*year,
#            family = poisson(link = "log"),data = deaths9311,
#            offset = log(ERP))

train=deaths9311 %>% slice_sample(prop=0.9)
test=anti_join(deaths9311,train)

# mod0081=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + Reg*year,
#             family = poisson(link = "log"),data = train,
#             offset = log(ERP))
rmse(mod001,mod002,mod003,mod004,mod005,mod006,mod007,mod008, mod009)
HelpersMG::compare_BIC(m1=mod001, m2=mod002,m3=mod003, m4=mod004, m5=mod005, m6=mod006, m7=mod007, m8=mod008,m9=mod009)
cv_mort09=cv_rmse(mod001,mod002,mod003,mod004,mod005,mod006,mod007,mod008,prop = 0.9,n = 10,data=deaths9311)



#fertility model selection ####
births_8116 = births_8116 %>%
  mutate(year=as_factor(year)) %>% 
  rename(x=births)

mod101=glm(round2(x,0)~(A+Reg+year),
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))
mod102=glm(round2(x,0)~A+Reg+year + A*year,
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))
mod103=glm(round2(x,0)~A+Reg+year + A*Reg + A*year,
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))
mod104=glm(round2(x,0)~A+Reg+year + A*Reg + A*year + Reg*year,
           family = poisson(link = "log"),data = births_8116,
           offset = log(ERP))


HelpersMG::compare_BIC(m1=mod101, m2=mod102, m3=mod103, m4=mod104)
rmse(mod101, mod102, mod103, mod104)


# immigration model selection ####
intmigr_8111 <- intmigr_8116 %>%
  mutate(x=immig,
         def=as_factor(year>2002) %>% fct_recode("old"="FALSE","new"="TRUE"),
         year=as.factor(year)
         ) 

#loglin models
mod201=glm(round2(x,0)~(A+S+Reg+year),
           family = poisson(link = "log"),data = intmigr_8111)
#AT
mod202=glm(round2(x,0)~A+S+Reg+year + A*year,
           family = poisson(link = "log"),data = intmigr_8111)
#AS 
mod203=glm(round2(x,0)~A+S+Reg+year + A*year + A*S,
           family = poisson(link = "log"),data = intmigr_8111)
#R effects
mod204=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg,
           family = poisson(link = "log"),data = intmigr_8111)
mod205=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + S*Reg,
           family = poisson(link = "log"),data = intmigr_8111)
mod206=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg + S*Reg,
           family = poisson(link = "log"),data = intmigr_8111)
#rmse() mod004
mod207=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg,
           family = poisson(link = "log"),data = intmigr_8111)
#S
# mod208=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg + S*year,
           # family = poisson(link = "log"),data = intmigr_8111)
mod208=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + S*Reg,
           family = poisson(link = "log"),data = intmigr_8111)

#Reg*year
# mod209=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg + Reg*year,
           # family = poisson(link = "log"),data = intmigr_8111)
mod209=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + Reg*year,
           family = poisson(link = "log"),data = intmigr_8111)

mod210=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + Reg*year+ S*Reg,
           family = poisson(link = "log"),data = intmigr_8111)

mod2041=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg +def,
           family = poisson(link = "log"),data = intmigr_8111)
mod2071=glm(round2(x,0)~A+S+Reg + A*S*year + A*Reg +def,
           family = poisson(link = "log"),data = intmigr_8111)
mod2042=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg +def*Reg,
            family = poisson(link = "log"),data = intmigr_8111)
mod2072=glm(round2(x,0)~A+S+Reg + A*S*year + A*Reg +def*Reg,
            family = poisson(link = "log"),data = intmigr_8111)

rmse(mod201,mod202,mod203,mod204,mod205,mod206,mod207,mod208,mod209,mod210)
HelpersMG::compare_BIC(m1=mod201, m2=mod202,m3=mod203, m4=mod204, m5=mod205, m6=mod206, m7=mod207, m8=mod208, m9=mod209, m10=mod210)
cv_immi09=cv_rmse(mod201,mod202,mod203,mod204,mod205,mod206,mod207,mod208, mod209,mod210,prop = 0.9,n = 10,data=intmigr_8116)

rmse(mod207,mod208,mod209,mod210,mod2041,mod2042,mod2071)

# for (i in 1:length(mods_mort)){
#   mod=glm(formula = mods_mort[[i]]$formula,
#           family = poisson(link = "log"),data = train,
#           offset = log(ERP))
#   preds=predict.glm(object = mod, newdata = test,type = "response")
#   print(sqrt(sum((preds-test$x)^2)))
# }




# emigration model selection ####
intmigr_8116 <- intmigr_8116 %>%
  mutate(x=emig,
         year=as.factor(year))


mod301=glm(round2(x,0)~(A+S+Reg+year),
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
#AT
mod302=glm(round2(x,0)~A+S+Reg+year + A*year,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
#AS 
mod303=glm(round2(x,0)~A+S+Reg+year + A*year + A*S,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
#R effects
mod304=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
mod305=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + S*Reg,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
mod306=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg + S*Reg,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
#rmse() mod004
mod307=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
#2way +ST
# mod308=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg + S*year,
#            family = poisson(link = "log"),data = intmigr_8116,
#            offset = log(ERP))
mod308=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + S*Reg,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))

mod309=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + Reg*year,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))

#2-way +RT
mod310=glm(round2(x,0)~A+S+Reg+year + A*year + A*S + A*Reg + Reg*year,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))
mod310=glm(round2(x,0)~A+S+Reg+year + A*S*year + A*Reg + Reg*year + S*Reg,
           family = poisson(link = "log"),data = intmigr_8116,
           offset = log(ERP))

rmse(mod301,mod302,mod303,mod304,mod305,mod306,mod307,mod308,mod309,mod310)
HelpersMG::compare_BIC(m1=mod301, m2=mod302,m3=mod303, m4=mod304, m5=mod305, m6=mod306, m7=mod307, m8=mod308, m9=mod309)


