# goodness of fit measures ####
#A unified framework for probabilistic forecasting subnational populations
#(C) anonymised (2020) 

# libraries ####
library(ggplot2)

# functions ####
# function fit.summary creates a df for plotting "residuals" and error bars against observed rates/data,
# it also produces a set of error GOF measures (ME, MAE, MAPE, RMSE)
fit.summary <- function(fit,data,co){
  dat=data %>% 
    filter(year<2012) %>% 
    {if (co == "f") select(.,births,ERP, Reg,A) %>% 
        mutate(xr=births/ERP) else .} %>%
    {if (co == "m") select(.,x,ERP, Reg,S,A) %>% 
        mutate(xr=x/ERP) else .} %>%
    {if (co == "im") select(.,immig, Reg,S,A) %>% 
        mutate(xr=immig) else .} %>%
    {if (co == "em") select(.,emig,ERP, Reg,S,A) %>% 
        mutate(xr=emig/ERP) else .} %>%
    {if (co == "in") filter(.,year<2011) %>% 
        select(Flow,x, Corridor,S,A) %>% 
        mutate(xr=Flow/x) else . 
    }  
  
  df=summary(fit,pars="mmd")$summary %>% 
    as.data.frame() %>% 
    bind_cols(dat) %>%
    mutate(dat_l=log(xr),
           resid=dat_l-`50%`)
  dff=df %>% filter(is.finite(resid))
  me=mean(dff$resid)
  mae=mean(abs(dff$resid))
  rmse=sqrt(mean(dff$resid^2))
  mape=dff %>% mutate(mape=abs(resid/dat_l)) %>% pull(mape) %>% mean()
  print(me)
  print(mae)
  print(mape)
  print(rmse)
  return(df)
}  

#calculation of the measures ####
sum_od084=fit.summary(fit08_04,input_au,co="in")
sum_im022=fit.summary(fit_imi022,intmigr_8116,co="im")
sum_im031=fit.summary(fit_imi031,intmigr_8116,co="im")
sum_em02=fit.summary(fit_emi02,intmigr_8116,co="em")
sum_f01=fit.summary(fit_f01,births_8116,co="f")
sum_f03=fit.summary(fit_f03,births_8116,co="f")
sum_m02=fit.summary(fit_m02,deaths9311,co="m")
sum_m03=fit.summary(fit_m03,deaths9311,co="m")

#total fit
plot.gof <- function(input,digits=0,save=T,low,high){
p=ggplot(data = input %>% rename(Sex=S)) + 
  geom_errorbar(aes(x=dat_l,ymin=`2.5%`,ymax=`97.5%`, colour=Sex), alpha=0.6) + 
  # geom_point(aes(x=immig_l,y=`50%`, colour=S),size=0.1) +
  geom_abline(intercept=0, slope=1) +
  theme_bw() +
  scale_x_continuous(limits=c(low,high),labels = function(x) format(round(exp(x),digits = digits), 
                                                 trim=T,digits=2)) + 
  scale_y_continuous(limits=c(low,high),labels = function(x) format(round(exp(x),digits = digits), 
                                                 trim=T,digits=2)) +
  # facet_wrap(A~.,ncol = 3, scales = "free") +
  labs(y="estimate (95% predictive interval)", x="data")
if (save==F) return(p)
if (save==T) {ggsave(filename = paste0("plots/Figure_fit_",deparse(substitute(input)),".pdf"), 
                    plot = p,
                    device="pdf",width = 3, 
                    height = 2, scale=2)}
}
range(sum_od084$`97.5%`)
range(sum_od084$`2.5%`)
range(sum_od084$dat_l)
plot.gof(sum_od084,2,low = -10,high=-1.5)
plot.gof(sum_im031,low = -5,high=10)
plot.gof(sum_em02,3,low = -12,high=-2.5)
sum_f01 %>% bind_cols(data.frame(S="F")) %>% plot.gof(digits = 5,low = -12,high=-1.5)
#plot.gof(sum_f03 %>% bind_cols(data.frame(S="F")),2)
plot.gof(sum_m02,3,low = -10,high=-1)


#fit by age
ggplot(data = sum_im031) + 
  geom_errorbar(aes(x=dat_l,ymin=`2.5%`,ymax=`97.5%`, colour=S)) + 
  # geom_point(aes(x=immig_l,y=`50%`, colour=S),size=0.1) +
  geom_abline(intercept=0, slope=1) +
  theme_bw() +
  scale_x_continuous(labels = function(x) format(round(exp(x),0), 
                                                 trim=T,digits=2)) + 
  scale_y_continuous(labels = function(x) format(round(exp(x),0), 
                                                 trim=T,digits=2)) +
  facet_wrap(A~.,ncol = 3, scales = "free") +
  labs(y="estimate (95% predictive interval)", x="data")



