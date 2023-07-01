###########################
# testing lasso approach to variable selection
# CONCLUSION: the most complex model with 2-way interactions is preferred which is overly parameterised approach with lack of variability in some Corridor*Age groups
# LASSO ####
## glmnet ####
library(glmnet)
x=data.matrix(input_au[,c('Reg_O','Reg_D','Corridor','A','S','year')])
x=model.matrix(mod3)
y=round2(input_au$Flow,0)
modlass3=glmnet(x, y, family = poisson(link = "log"),alpha=1,offset = log(input_au$x))
coef(modlass2)
cv_model <- cv.glmnet(x, y, family = poisson(link = "log"),alpha=1,offset = log(input_au$x))
best_lambda <- cv_model$lambda.min
best_lambda
plot(cv_model)
modlass2=glmnet(x, y, family = poisson(link = "log"),alpha=1,offset = log(input_au$x), lambda = best_lambda)
coef(modlass2)
#mod2=glm(round2(Flow,0)~Corridor+year+A*Reg_O+A*Reg_D+S*A-1,family = poisson(link = "log"),data = input_au,offset = log(x))

#######################
## gglasso ####
library(gglasso)
input_au1 = input_au %>%
  mutate(year=as.factor(year))
y=log(((input_au1$Flow+1)*100000)/input_au1$x)
x1 <- model.matrix( ~ (A+S+year+Corridor)^2, dplyr::select(input_au1, A,S,year,Corridor))[,-1]
colnames(x1)
groups=c(rep(1,16),2,rep(3,6),rep(4:10,each=16),rep(11,6))
#with intercept
groups=c(0,rep(1,16),2,rep(3,6),rep(4,55),rep(5,16), rep(6,16*6),
         rep(7,16*55), rep(8,6), rep(9,55), rep(10,55*6))+1
#without intercept
groups=c(rep(1,16),2,rep(3,6),rep(4,55),rep(5,16), rep(6,16*6),
         rep(7,16*55), rep(8,6), rep(9,55), rep(10,55*6))
fit1 <- gglasso(x = x1, y = y, group = groups, lambda = 1)
fit0027 <- gglasso(x = x1, y = y, group = groups, lambda = 0.02779054)
#fit$beta
fit0027$beta %>% as.data.frame() %>% filter(!s0==0) %>% View
fit_cv <- gglasso(x = x1, y = y, group = groups, nlambda = 30)

fit_cv22 <- cv.gglasso(x = x1, y = y, group = groups, lambda = NULL,
                       pred.loss = "L1", nfolds = 5)
fit_cv2 <- cv.gglasso(x = x1, y = y, group = groups, lambda = NULL,
                      pred.loss = "L1", nfolds = 5)
#using log-rates
fit11_0 <- gglasso(x = x1, y = y, group = groups, lambda = fit_cv11$lambda.min)
fit11_0$beta %>% as.data.frame() %>% filter(!s0==0) %>% View
