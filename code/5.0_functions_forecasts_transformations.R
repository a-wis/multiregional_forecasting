# transofmrations of raw results into data frames for easy plotting
#A unified framework for probabilistic forecasting subnational populations
#(C) anonymised (2020) 

# generating forecasts ####
#this function generates forecasts based on the output from stan - it requires sampling from the normal distribution - this is way much faster than doing it in stan
gen.forecast=function(inp=fit_f0001,S=1,R=8,Fh=15,Ti=31,Ag=7,start.c=1,fit.too=T)
{
  #inp = stanfit
  require(rstan)
  fit.matrix=as.matrix(inp)
  size=dim(fit.matrix)
  #reading in std of log-rates
  sig1 = fit.matrix[start.c:size[1],which(colnames(fit.matrix)=="sig1")]
  #reading in expected log-rates
  rds=NULL
  if (fit.too==T){
    fit.rate=fit.matrix[start.c:size[1],which(colnames(fit.matrix)=="mmd[1]"):which(colnames(fit.matrix)==paste0("mmd[",R*Ag*S*Ti,"]"))]
    rds=t(matrix(rnorm((size[1]-start.c+1)*R*Ag*(Ti)*S,fit.rate,sig1),size[1]-start.c+1,R*Ag*(Ti)*S))
  }
  #reading in forecasts of log-rates
  forc.rate = fit.matrix[start.c:size[1],which(colnames(fit.matrix)=="mmd_f[1]"):which(colnames(fit.matrix)==paste0("mmd_f[",R*Ag*S*Fh,"]"))]  #AxT
  #sampling with std sig1
  rfs=t(matrix(rnorm((size[1]-start.c+1)*R*Ag*(Fh)*S,forc.rate,sig1), size[1]-start.c+1,R*Ag*(Fh)*S))
  #binding fit with forecasts 
  return(rbind(rds,rfs)) #rfit=fit.rate.f,#abind(fit.rate.f,rf,along=3)
}

#using mpd parameter for fit UNUSED ####
gen.forecast.p=function(inp=fit_f0001,dat,S=1,R=8,Fh=15,Ti=31,Ag=7,start.c=1,fit.too=T)
{
  #inp = stanfit
  require(rstan)
  fit.matrix=as.matrix(inp)
  size=dim(fit.matrix)
  #reading in std of log-rates
  sig1 = fit.matrix[start.c:size[1],which(colnames(fit.matrix)=="sig1")]
  #reading in expected log-rates
  rds=NULL
  if (fit.too==T){
    fit.rate=t(fit.matrix[start.c:size[1],which(colnames(fit.matrix)=="mpd[1]"):which(colnames(fit.matrix)==paste0("mpd[",R*Ag*S*Ti,"]"))])
    rds=fit.rate-log(dat)
  }
  #reading in forecasts of log-rates
  forc.rate = fit.matrix[start.c:size[1],which(colnames(fit.matrix)=="mmd_f[1]"):which(colnames(fit.matrix)==paste0("mmd_f[",R*Ag*S*Fh,"]"))]  #AxT
  #sampling with std sig1
  rfs=t(matrix(rnorm((size[1]-start.c+1)*R*Ag*(Fh)*S,forc.rate,sig1), size[1]-start.c+1,R*Ag*(Fh)*S))
  #binding fit with forecasts 
  return(rbind(rds,rfs)) #rfit=fit.rate.f,#abind(fit.rate.f,rf,along=3)
}



#transforming result into data frame ####
transform.df<-function(mat_obj=f16_mat,dat=births_8116,yd0=1981,ydl=2011,yf0=2012,yfl=2026,pop.comp="F", interval=5){
  require(dplyr)
  if (pop.comp=="F"){
    temp0a = list(Reg=1:8,year=yd0:ydl,A=levels(dat$A))
    temp0b = list(Reg=1:8,year=yf0:yfl,A=levels(dat$A))
  } 
  if (pop.comp=="OD"){
    temp0a = list(Reg_D=1:8,Reg_O=1:8,year=seq(yd0,ydl,interval),S=levels(dat$S),A=levels(dat$A)) # yd0=1981, ydl=2006
    temp0b = list(Reg_D=1:8,Reg_O=1:8,year=seq(yf0+4,yfl,interval),S=levels(dat$S),A=levels(dat$A)) # yf0=2007, yfl=2021
  }
  if (pop.comp%in%c("M","I","E")){
    temp0a = list(Reg=1:8,year=yd0:ydl,S=levels(dat$S),A=levels(dat$A))
    temp0b = list(Reg=1:8,year=yf0:yfl,S=levels(dat$S),A=levels(dat$A))
  }
  temp1a=expand.grid(temp0a)
  temp1b=expand.grid(temp0b)
  temp1=bind_rows(temp1a,temp1b)
  if (pop.comp!="OD") temp1=temp1 %>% mutate(Reg=as.factor(Reg))
  if (pop.comp=="OD") temp1=temp1 %>% mutate(Corridor = paste0(Reg_O,"->",Reg_D)) %>%
    filter(Reg_O != Reg_D, A != 0) 
  result=cbind(temp1,exp(mat_obj))
  result=gather(data = result, "Iter", "x","1":"1000")
  return(result)
}


##### life tables ####

#creating life table for estimates ####
#mortality rate
life.tab.f<-function(inp1, Ti=15){
  #if (gfit==F){
  inm=aperm(array(t(exp(inp1)),c(1000,8,Ti,2,18)),c(1,2,5,4,3))#rast
  # inm=inm[,,,,c(3,8,13)]
  #   L=30
  #   } else {
  #     inm=array(inp1[,which(colnames(inp1)=="md[1,1,1,1]"):which(colnames(inp1)=="md[5,16,2,5]")],c(500,5,16,2,5))#rast
  #     L=5
  #   }
  I=array(0,c(1000,8,8,18,2,Ti)) 
  I=with(expand.grid(a = 1:1000, b = 1:8,c=1:18,d=1:2,w=1:Ti), replace(I, cbind(a, b, b,c,d,w), 1))
  M1=array(0,c(1000,8,8,18,2,Ti))
  
  M1=with(expand.grid(a = 1:1000, b = 1:8,c=1:18,d=1:2,w=1:Ti), replace(M1, cbind(a, b, b,c,d,w), inm))#+apply(od1,c(1:2,4:6),sum)))-od1#+apply(od,c(1:2,4:6),sum)-od
  s2=aperm(array(apply(I+2.5*M1,c(1,4,5,6),solve),c(8,8,1000,18,2,Ti)),c(3,1,2,4,5,6))
  
  s21=I-2.5*M1
  #survival prob
  
  #here equation 22 is used
  p=array(NA,c(1000,8,8,18,2,Ti))
  
  for (i in 1:1000){
    for (j in 1:2){
      for (t in 1:Ti){
        for (a in 1:18){
          p[i,,,a,j,t]=s2[i,,,a,j,t]%*%s21[i,,,a,j,t]}
      } } }
  
  
  #Here Rogers 1975 is used, page 65
  l=array(0,c(1000,8,8,18,2,Ti)) #persons
  L=array(0,c(1000,8,8,18,2,Ti)) #person-years
  lradix= I*100
  #radix
  for (i in 1:1000){
    for (j in 1:2){
      for (t in 1:Ti){
        l[i,,,1,j,t]=p[i,,,1,j,t]%*%lradix[i,,,a,j,t]
        for (a in 1:17){ 
          l[i,,,a+1,j,t]=p[i,,,a,j,t]%*%l[i,,,a,j,t]
          L[i,,,a,j,t]=2.5*(I[i,,,a,j,t]+p[i,,,a,j,t])%*%l[i,,,a,j,t]
        }
        L[i,,,18,j,t]=solve(M1[i,,,18,j,t])%*%l[i,,,18,j,t]
      } } }
  
  #person years left to live
  Tt=array(0,c(1000,8,8,18,2,Ti))
  Tt[,,,18,,]=L[,,,18,,]
  
  for (a in 1:17){
    Tt[,,,a,,]=apply(L[,,,a:18,,],c(1,2,3,5,6),sum)
  }
  #life expectancies
  lex=array(0,c(1000,8,8,18,2,Ti))
  
  for (i in 1:1000){
    for (j in 1:2){
      for (t in 1:Ti){
        for (a in 1:18){ 
          lex[i,,,a,j,t]=Tt[i,,,a,j,t]%*%solve(l[i,,,a,j,t])
        }}}}
  return(lex)
  
}




#life table for data ####
#mortality rate
life.tab.d<-function(){
  md=aperm(array(deaths9311$x/deaths9311$ERP,c(8,19,2,18)),c(1,4,3,2))
  I=array(0,c(1000,8,8,18,2,19)) 
  I=with(expand.grid(a = 1:1000, b = 1:8,c=1:18,d=1:2,w=1:19), replace(I, cbind(a, b, b,c,d,w), 1))
  Md1=array(0,c(8,8,18,2,19))
  Md1=with(expand.grid(b = 1:8,c=1:18,d=1:2,w=1:19), replace(Md1, cbind(b, b,c,d,w), md))#+apply(odd[,,,,13:17],c(1,3:5),sum)))-odd[,,,,13:17]#+apply(odd[,,,,13:17],c(1,3:5),sum)-odd[,,,,13:17]
  s2=array(apply(I[1,,,,,1:19]+2.5*Md1,c(3,4,5),solve),c(8,8,18,2,19))
  s21=I[1,,,,,1:19]-2.5*Md1
  #survival prob
  
  #here equation 22 is used
  pd=array(NA,c(8,8,18,2,19))
  
  for (j in 1:2){
    for (t in 1:19){
      for (a in 1:18){
        pd[,,a,j,t]=s2[,,a,j,t]%*%s21[,,a,j,t]}
    }}
  
  
  #Here Rogers 1975 is used, page 65
  ld=array(0,c(8,8,18,2,19)) #persons
  Ld=array(0,c(8,8,18,2,19)) #person-years
  ldradix= I[1,,,,,1:19]*100
  #radix
  for (j in 1:2){
    for (t in 1:19){
      ld[,,1,j,t]=pd[,,1,j,t]%*%ldradix[,,a,j,t]
      for (a in 1:17){ 
        ld[,,a+1,j,t]=pd[,,a,j,t]%*%ld[,,a,j,t]
        Ld[,,a,j,t]=2.5*(I[1,,,a,j,t]+pd[,,a,j,t])%*%ld[,,a,j,t]}
      Ld[,,18,j,t]=solve(Md1[,,18,j,t])%*%ld[,,18,j,t]
    }}
  
  #person years left to live
  Ttd=array(0,c(8,8,18,2,19))
  Ttd[,,18,,]=Ld[,,18,,]
  
  for (a in 1:17){
    Ttd[,,a,,]=apply(Ld[,,a:18,,],c(1,2,4,5),sum)
  }
  #life expectancies
  lexd=array(0,c(8,8,18,2,19))
  
  
  for (j in 1:2){
    for (t in 1:19){
      for (a in 1:18){ 
        lexd[,,a,j,t]=Ttd[,,a,j,t]%*%solve(ld[,,a,j,t])
      }}}
  
  return(lexd)
}

# mortality results transfomation ####
#transforms the mortality rates into outputs from the life table - life expectancies
transform.mort.lex<-function(mat_m){
  mat_m.df=transform.df(mat_obj = mat_m,dat = deaths9311,yd0 = 1993,ydl = 2011,pop.comp = "M")
  mat_m.df=arrange(mat_m.df,A,S,year,Reg)
  t_m04=life.tab.f(log(mat_m.df$x),34)
  t_m04=aperm(apply(t_m04,c(1,4,5,6),diag),c(2,1,3,4,5))
  t_m04=as.data.frame.table(t_m04)
  t_m04 = t_m04 %>% rename(Iter=Var1, Reg=Var2, A=Var3, S=Var4, year=Var5,x=Freq) %>%
    mutate(Iter=as.character(plyr::mapvalues(Iter,levels(Iter),1:1000)),
           Reg=plyr::mapvalues(Reg,levels(Reg),Regn),
           A=plyr::mapvalues(A,levels(A),levels(deaths9311$A)),
           S=plyr::mapvalues(S,levels(S),levels(deaths9311$S)),
           year=as.integer(1992+as.integer(plyr::mapvalues(year,levels(year),c(1993:2026)))))
  return(t_m04)
}