#making interregional projections
#A unified framework for probabilistic forecasting subnational populations
#(C) anonnymised (2020) 

#mortality, int emi- and immigration rates need to be arrays of size
#(#iterations, #regions, #age-groups, #sexes, #years)
#all should be in levels (not logs)

#fertility rates size is #(#iterations, #regions, #age-groups, #years)
#levels, not logs

#internal migration should be #(iterations, destination, origin, age-1, sex, time)
#levels, not logs
#with ZEROS on origin-destination DIAGONAL

  ##make a migration survivorship matrix S_bar (4.25) 
  # last age group is half of second to last
int.migr.s=function(inp=inp_od1){
    fun1=function(x){diag(x)<-1-colSums(x); return(x)}  
    temp=aperm(array(apply(inp,c(1,4,5,6),fun1),c(8,8,1000,18,2,3)),c(3,1,2,4,5,6))
    return(temp)
  }
  #P_bar (4.25) first part


# comp.od() function transforms internal migration results into a transition matrix
comp.od<-function(inp=mat_od){
  temp=aperm(array(t(exp(inp)),c(1000,7,8,3,2,17)),c(1,2,3,6,5,4))#doast
  od=array(0,c(1000,8,8,17,2,3))
  for (j in 1:7){
    for (i in 1:j){ #upper triangle
      od[,i,j+1,,,]=temp[,i,j+1,,,]}
    for (i in (j+1):8){ # this is lower triangle
      od[,i,j,,,]=temp[,i-1,j,,,]}
  }
  od1=od#aperm(od,c(1,3,2,4,5,6)) # IT IS ALREADY TRANSPOSED
  od2=array(0,c(1000,8,8,18,2,3))
  od2[,,,1:17,,]=od1 #that's what we need
  #array(apply(od,c(1,4,5,6),t),c(8,8,1000,18,2,3)),c(3,1,2,4,5,6))
  return(od2)
}

 # pop.for.op() function forecasts pop using inputs on mortality, fertility, immigration, emigration, internatl migration, with specfied Forecast horison Fh.i (default three periods). It uses 5-year age groups and offers little to no flexibility about other age groups. 
pop.for.op<-function(m,f,im,em,od1,Fh.i=3){ 
  #m - mortality rates d (in manuscript), f- fertility rates, im - immigration counts (manuscript G),
  #em - emigration rates (manusctript e), od1 - internal migration probabilities (manuscript m)
#computing survivorship s
    It=dim(m)[1]
    Re=dim(m)[2]
    Ag=dim(m)[3]
    Se=dim(m)[4]
    Fh=dim(m)[5]
    Fhf=Fh.i
    # M=array(0,c(500,5,5,16,2,30))#it,reg,reg,ag,sex,time
    # S_bar (manuscript Eq.5)
    S_bar=array(0,c(It,Re,Re,Ag,Se,Fh.i))#it,reg,reg,ag,sex,time
    S_bar=int.migr.s(od1) #(4.25)
    #(manuscript Eq.4)
    P_bar=array(0,c(It,Re,Re,Ag,Se,Fh.i))#it,reg,reg,ag,sex,time
    P_bar[,,,2:(Ag-1),,] = .5*S_bar[,,,1:(Ag-2),,]+.5*S_bar[,,,2:(Ag-1),,]
    P_bar[,,,Ag,,] = S_bar[,,,Ag-1,,] #last age group
    P_bar[,,,1,,] = .5*aperm(array(apply(S_bar[,,,1,,],c(1,4,5),function(x){x%*%x}),c(Re,Re,It,2,Fh.i)),c(3,1,2,4,5)) + 0.5*S_bar[,,,1,,] #4.28
    
    #(manuscript Eq.4)
    P_delta=array(0,c(It,Re,Re,Ag,Se,Fh.i))#it,reg,reg,ag,sex,time #(manuscript Eq.6)
    P_bar_sum=array(0,c(It,Re,Re,Ag,Se,Fh))
    for (i in 1:It){
      for (j in 1:2){
        for (t in 1:Fh){
          for (a in 1:Ag){
            P_bar_sum[i,,,a,j,t] = P_bar[i,,,a,j,(t+4)%/%5]*(m[i,,a,j,t]+em[i,,a,j,t])
          }
        }
      }
    }
    P_bar_sum1=unname(apply(P_bar_sum,c(1:5),tapply, c(rep(1:Fh.i,each=5)), sum))
    P_bar_sum2=aperm(apply(P_bar_sum1,c(1,2,5,6),function(x){1+0.5*colSums(x)}),c(2,3,1,4,5))
    m_em1=1-0.5*unname(apply(m+em,c(1:4),tapply, c(rep(1:Fh.i,each=5)), sum))
    P_bar_sum3=m_em1/P_bar_sum2
    for (i in 1:It){
      for (j in 1:2){
        for (t in 1:Fh.i){
          for (a in 1:Ag){
            for (r in 1:Re){
            P_delta[i,r,r,a,j,t]=P_bar_sum3[t,i,r,a,j] 
    }}}}}
    ##
    #M from 4.31 (Rogers 1995)
    # M=array(0,c(It,Re,Re,Ag,S,Fh.i))#it,reg,reg,ag,sex,time
    I=array(0,c(It,Re,Re,Ag,Se,Fh.i))
    I=with(expand.grid(a = 1:It, b = 1:Re,c=1:Ag,d=1:2,w=1:Fh.i), replace(I, cbind(a, b, b,c,d,w), 1))
    
    # Fr is fertility B in manuscript Eq.1
    Fr=array(0,c(It,Re,Re,9,Fh))
    Fr=with(expand.grid(a = 1:It, b = 1:Re,c=1:9,w=1:Fh), replace(Fr, cbind(a, b, b,c,w), f))
 print("Now calculating matrix P and S")   
    ##
 S=array(0,c(It,Re,Re,Ag,2,Fh.i)) #matrix S=(I+P(x)) %*% P(x) %*% (I+P(x))^-1; p.101 Rogers 1995
 P=array(0,c(It,Re,Re,Ag,2,Fh.i)) #matrix P=P_bar %*% P_delta; p.100-101 (4.26 Rogers 1995) or manuscript Eq.4
 S_n=array(0,c(It,Re,Re,2,Fh.i)) #matrix survivorship of infants (manuscript Eq.7)
 for (i in 1:It){
   for (j in 1:2){
     for (t in 1:Fh.i){
       for (a in 1:Ag){
         P[i,,,a,j,t]=P_bar[i,,,a,j,t]%*%P_delta[i,,,a,j,t]
         # M[i,,,a,j,t]=.4*solve(I[i,,,a,j,t]+P[i,,,a,j,t])%*%(I[i,,,a,j,t]-P[i,,,a,j,t])
       }
       for (a in 1:(Ag-1)){
        S[i,,,a,j,t]=(I[i,,,a,j,t]+P[i,,,a+1,j,t])%*%P[i,,,a,j,t]%*%solve(I[i,,,a,j,t]+P[i,,,a,j,t])
       }
       S[i,,,Ag,j,t]=P[i,,,Ag,j,t]
         # I[i,,,Ag,j,t]-(I[i,,,Ag,j,t]+P[i,,,Ag,j,t])%*%solve(I[i,,,Ag,j,t]+P[i,,,Ag,j,t])%*%(I[i,,,Ag,j,t]-P[i,,,Ag,j,t]) # last group surv Option 1 (4.35)
       S_n[i,,,j,t]=solve(I[i,,,1,j,t]+solve(I[i,,,1,j,t]+P[i,,,1,j,t])%*%(I[i,,,1,j,t]-P[i,,,1,j,t])) #2nd Eq. p 101, then the approximation for newborns S(-5)
     }}}
 
 #population baseline - 
 K=array(0,c(Re,Ag+1,2)) #Ag+1 to make it easier to compute newborns
 K[,1:Ag,]=aperm(array(filter(ERP_pop,year==2011)$ERP,c(Ag,2,Re)),c(3,1,2))
 K[,Ag+1,]=K[,Ag,]
 
 print("... first year of forecast") 
  #Pf- forecasted pop
 #manuscript Eqs. 1 & 2 for births i.e. extended with international migration
 #first year (5-year interval)
 Pf=array(NA,c(It,Re,Ag,2,Fh.i))
 Fra=array(NA,c(It,Re,Fh.i))
 Fr1=array(NA,c(It,Re,2,Fh.i))
 for (i in 1:It){
   for (t in 1:1){
     for (j in 1:2){
       for (a in 1:(Ag-2)){
         Pf[i,,a+1,j,t]=S[i,,,a,j,t]%*%(K[,a,j]+apply(im[i,,a,j,((5*t-4):(5*t))],c(1),sum)*.5)+apply(im[i,,a+1,j,((5*t-4):(5*t))],c(1),sum)*.5
         }
      Pf[i,,Ag,j,t]=S[i,,,Ag-1,j,t]%*%(K[,Ag-1,j] + apply(im[i,,Ag-1,j,((5*t-4):(5*t))],c(1),sum)*.25+apply(im[i,,Ag,j,((5*t-4):(5*t))],c(1),sum)*.25) + S[i,,,Ag,j,t]%*%K[,Ag,j] + apply(im[i,,Ag,j,((5*t-4):(5*t))],c(1),sum)*.5
      }
  Fra[i,,t]=(Fr[i,,,1,5*t]+Fr[i,,,2,5*t]%*%S[i,,,3,2,t])%*%(K[,3,2]+apply(im[i,,3,2,((5*t-4):(5*t))],c(1),sum)*.5)+
  (Fr[i,,,2,5*t]+Fr[i,,,3,5*t]%*%S[i,,,4,2,t])%*%(K[,4,2]+apply(im[i,,4,2,((5*t-4):(5*t))],c(1),sum)*.5)+
  (Fr[i,,,3,5*t]+Fr[i,,,4,5*t]%*%S[i,,,5,2,t])%*%(K[,5,2]+apply(im[i,,5,2,((5*t-4):(5*t))],c(1),sum)*.5)+
  (Fr[i,,,4,5*t]+Fr[i,,,5,5*t]%*%S[i,,,6,2,t])%*%(K[,6,2]+apply(im[i,,6,2,((5*t-4):(5*t))],c(1),sum)*.5)+ 
  (Fr[i,,,5,5*t]+Fr[i,,,6,5*t]%*%S[i,,,7,2,t])%*%(K[,7,2]+apply(im[i,,7,2,((5*t-4):(5*t))],c(1),sum)*.5)+ 
  (Fr[i,,,6,5*t]+Fr[i,,,7,5*t]%*%S[i,,,8,2,t])%*%(K[,8,2]+apply(im[i,,8,2,((5*t-4):(5*t))],c(1),sum)*.5)+
  (Fr[i,,,7,5*t]+Fr[i,,,8,5*t]%*%S[i,,,9,2,t])%*%(K[,9,2]+apply(im[i,,9,2,((5*t-4):(5*t))],c(1),sum)*.5)+
(Fr[i,,,8,5*t]+Fr[i,,,9,5*t]%*%S[i,,,10,2,t])%*%(K[,10,2]+apply(im[i,,10,2,((5*t-4):(5*t))],c(1),sum)*.5)
  Fr1[i,,1,t]=0.5121951*Fra[i,,t]#sexes split
  Fr1[i,,2,t]=0.4878049*Fra[i,,t]
         Pf[i,,1,1,t]=2.5*S_n[i,,,1,t]%*%(Fr1[i,,1,t])+apply(im[i,,1,1,((5*t-4):(5*t))],c(1),sum)*.25
         Pf[i,,1,2,t]=2.5*S_n[i,,,2,t]%*%(Fr1[i,,2,t])+apply(im[i,,1,2,((5*t-4):(5*t))],c(1),sum)*.25
     }}
 print("... and the rest of years")   
 #the rest of the years
 for (i in 1:It){
   for (t in 2:Fh.i){
     for (j in 1:2){
       for (a in 1:(Ag-2)){
         Pf[i,,a+1,j,t]=S[i,,,a,j,t]%*%(Pf[i,,a,j,t-1]+apply(im[i,,a,j,((5*t-4):(5*t))],c(1),sum)*.5)+apply(im[i,,a+1,j,((5*t-4):(5*t))],c(1),sum)*.5
       }
         Pf[i,,Ag,j,t]=S[i,,,Ag-1,j,t]%*%(Pf[i,,Ag-1,j,t-1]+apply(im[i,,Ag-1,j,((5*t-4):(5*t))],c(1),sum)*.25+apply(im[i,,Ag,j,((5*t-4):(5*t))],c(1),sum)*.25)+apply(im[i,,Ag,j,((5*t-4):(5*t))],c(1),sum)*.5 + S[i,,,Ag,j,t]%*%Pf[i,,Ag,j,t-1]}
     Fra[i,,t]=(Fr[i,,,1,5*t]+Fr[i,,,2,5*t]%*%S[i,,,3,2,t])%*%(Pf[i,,3,1,t-1]+apply(im[i,,3,2,((5*t-4):(5*t))],c(1),sum)*.5)+
(Fr[i,,,2,5*t]+Fr[i,,,3,5*t]%*%S[i,,,4,2,t])%*%(Pf[i,,4,1,t-1]+apply(im[i,,4,2,((5*t-4):(5*t))],c(1),sum)*.5)+
(Fr[i,,,3,5*t]+Fr[i,,,4,5*t]%*%S[i,,,5,2,t])%*%(Pf[i,,5,1,t-1]+apply(im[i,,5,2,((5*t-4):(5*t))],c(1),sum)*.5)+
(Fr[i,,,4,5*t]+Fr[i,,,5,5*t]%*%S[i,,,6,2,t])%*%(Pf[i,,6,1,t-1]+apply(im[i,,6,2,((5*t-4):(5*t))],c(1),sum)*.5)+ 
(Fr[i,,,5,5*t]+Fr[i,,,6,5*t]%*%S[i,,,7,2,t])%*%(Pf[i,,7,1,t-1]+apply(im[i,,7,2,((5*t-4):(5*t))],c(1),sum)*.5)+ 
(Fr[i,,,6,5*t]+Fr[i,,,7,5*t]%*%S[i,,,8,2,t])%*%(Pf[i,,8,1,t-1]+apply(im[i,,8,2,((5*t-4):(5*t))],c(1),sum)*.5)+
(Fr[i,,,7,5*t]+Fr[i,,,8,5*t]%*%S[i,,,9,2,t])%*%(Pf[i,,9,1,t-1]+apply(im[i,,9,2,((5*t-4):(5*t))],c(1),sum)*.5)+
(Fr[i,,,8,5*t]+Fr[i,,,9,5*t]%*%S[i,,,10,2,t])%*%(Pf[i,,10,1,t-1]+apply(im[i,,10,2,((5*t-4):(5*t))],c(1),sum)*.5)
     Fr1[i,,1,t]=0.5121951*Fra[i,,t]#sexes split
     Fr1[i,,2,t]=0.4878049*Fra[i,,t]     
     Pf[i,,1,1,t]=2.5*S_n[i,,,1,t]%*%(Fr1[i,,1,t])+apply(im[i,,1,1,((5*t-4):(5*t))],c(1),sum)*.25
     Pf[i,,1,2,t]=2.5*S_n[i,,,2,t]%*%(Fr1[i,,2,t])+apply(im[i,,1,2,((5*t-4):(5*t))],c(1),sum)*.25
     }}

 return(Pf)
  }
  
  
