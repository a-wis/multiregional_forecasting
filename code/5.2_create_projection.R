#making interregional population projections
#A unified framework for probabilistic forecasting subnational populations
#(C) anonnymised (2020) 

#reading in data ####
source("code/data_processing.R")
#reading in source files with functions ####
source("code/functions_forecasts_transformations.R")
source("code/functions_multiregional_projection.R")
# source("code/component_estimation.R") #this will take long time (about 3-11 hours per pop component)

#transforming forecasts of components ####
mat_m=gen.forecast(inp = fit_m02,S=2,R=8,Fh = 15,Ti=19,Ag=18,fit.too =F)
mat_f=gen.forecast(inp = fit_f01,S=1,R=8,Fh = 15,Ti=31,Ag=7,fit.too = F)
mat_im=gen.forecast(inp = fit_imi031,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =F)
mat_em=gen.forecast(inp = fit_emi03,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =F)
mat_od=gen.forecast(inp = fit08_04,S=2,R=56,Fh = 3,Ti=6,Ag=17,fit.too =F)

#old 2020
# mat_m=gen.forecast(inp = fit_m0124,S=2,R=8,Fh = 15,Ti=19,Ag=18,fit.too =F)
# mat_f=gen.forecast(inp = fit_f0016,S=1,R=8,Fh = 15,Ti=31,Ag=7,fit.too = F)
# mat_im=gen.forecast(inp = fit_imi0044,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =F)
# mat_em=gen.forecast(inp = fit_emi0044,S=2,R=8,Fh = 15,Ti=31,Ag=18,fit.too =F)
# mat_od=gen.forecast(inp = fit03111,S=2,R=56,Fh = 3,Ti=6,Ag=17,fit.too =F)
# changing dimensions ####
#mortality 
inp_m4<-aperm(array(t(exp(mat_m)),c(1000,8,15,2,18)),c(1,2,5,4,3))#rast
#fertility  
f1<-aperm(array(t(exp(mat_f)),c(1000,8,15,7)),c(1,2,4,3))#rat #[1737:2576,]
inp_f16=array(0,c(1000,8,9,15))
inp_f16[,,2:8,]=f1
rm(f1)
#international immigration 02
inp_im04=aperm(array(t(exp(mat_im)),c(1000,8,15,2,18)),c(1,2,5,4,3))#rast  #international emigration  
inp_em04=aperm(array(t(exp(mat_em)),c(1000,8,15,2,18)),c(1,2,5,4,3))#rast
#internal migration
inp_od1=comp.od(inp=mat_od)

# creating projections ####
require(stats)
system.time(Pf_test<-pop.for.op(inp_m4,inp_f16,inp_im04,inp_em04,inp_od1,Fh.i = 3))
