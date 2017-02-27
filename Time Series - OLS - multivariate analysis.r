# # Assignment 2: Model Selection
rm(list = ls())
cat("\f")

# install.packages("glmnet")
#Loading packages
library(rJava)
library(xlsx)
library(glmnet)

#Define OLS
OLS <- function(X0,y0,varargin){
  # t(x0) is transpose of X0
  Beta = solve(t(X0)%*%X0)%*%t(X0)%*%y0
  yhat = X0%*%Beta        # fitted y
  nobs = nrow(X0)    
  nreg = ncol(X0)
  res = y0-yhat                          #Residuals
  shat = sum(res^2)/(nobs-nreg)          #variance of residuals
  varcovar = shat*solve(t(X0)%*%X0)      #var-covar matrix
  tstat = Beta/ sqrt(diag(varcovar))
  ESS = sum(res^2)
  TSS = sum((y-mean(y))^2)
  Rq = 1-ESS/TSS                         #Rq
  stilda = ESS/(nobs-(nreg-1)-1)
  S2 = (TSS)/(nobs-1)
  Rqad = 1-stilda/S2                               #Rqadj (maximize)
  Aic = log((ESS/nobs))+(2*nreg)/nobs              #Akaike (minimize)
  Bic = log((ESS/nobs))+nreg*log(nobs)/nobs        #Bic    (minimize)
  HQ = log((ESS/nobs))+2*nreg*log(log(nobs))/nobs  #HQ     (minimize)
  
  S = list(Beta,yhat,tstat,res,shat,varcovar,Rq,Rqad,Aic,Bic,HQ)
  names(S) = c('Beta','yhat','tstat','res','varRes','varcovarBeta','Rq','Rqadj','Akaike','Bic','HQ')
  
  S
}


#Define Forecasts2K
Forecasts2K <- function(X0,y0,Xfor0){
  # This function generates all possible 2^K forecasts using a TxK matrix of
  # predictors, X, and a Tx1 univariate dependent variable, y
  # Note the timing: X(t-1,:), y(t), Xfor(t,:) go in each row
  # X does not contain an intercept, but an intercept is always included in the
  # regression
  
  # X0 <- X[1 : 119]
  # y0 <- Y[2 : 120]
  # Xfor0 <- X[120]
  Nvar = ncol(X0)  #number of X-variables
  # Say for 3 independent x variables - there are 8 possibilities of inclusion or exclusion of vars.
  # This is what is counted below.
  # AllModels
  # X1 X2 X3
  # [1,]  0  0  0
  # [2,]  1  0  0
  # [3,]  0  1  0
  # [4,]  1  1  0
  # [5,]  0  0  1
  # [6,]  1  0  1
  # [7,]  0  1  1
  # [8,]  1  1  1
  
  
  AllModels = as.matrix(expand.grid(data.frame(rbind(numeric(Nvar),numeric(Nvar)+1))))  # all permutations of X-variables
  Criteria = matrix(NaN,nrow(AllModels),2) 
  AllF = rep(NaN,nrow(AllModels))
  intercept = rep(1,length(y0)) # 119
  nmodels = nrow(AllModels)      #number of models; will return 8 for 3 independent variables
  
  for (j in 1:nmodels) {
    ## For all possible models find which independent var will have a positive coefficient in the model
    model = which((AllModels[j,])!=0)
    
    REG = OLS(cbind(intercept,X0[,model]),y0)
    Criteria[j,1] = REG$Akaike
    Criteria[j,2] = REG$Bic
    AllF[j] = t(c(1,Xfor0[model]))%*% REG$Beta
  }
  S = list(AllF,Criteria)
  names(S) = c('AllF','Criteria')
  S
}


#Readin data
datas = read.xlsx('Goyal_Welch_data.xlsx',1,colIndex = 1:20)
temp =sapply(datas, is.nan)
# head(temp)
datas = as.matrix(datas)
# [-(1:max(sum(temp))),]
nrow(datas)

X = datas[seq(1,nrow(datas)),c(6:15)] # use data across all years and 10 predictor variables, you need to change this
y = datas[seq(1,nrow(datas)),5]      #   you need to change this
estimationEnd = 516 #% amount of data to use for initial forecast, 43 years -> 1927 until 1969 here, you need to change this
TT = nrow(X) # sample size

#Initialize data structures

  modelLabels = c('Data','AIC','BIC','Forward stepwise','Backward stepwise','Lasso','Combination','Kitchen Sink','PM')
  
  AICforecast = matrix(0,TT-estimationEnd,1)
  BICforecast = matrix(0,TT-estimationEnd,1)
  forwardStepForecast = matrix(0,TT-estimationEnd,1)
  backStepForecast = matrix(0,TT-estimationEnd,1)
  lassoForecast = matrix(0,TT-estimationEnd,1)
  combinationForecast = matrix(0,TT-estimationEnd,1)
  kitchenSink = matrix(0,TT-estimationEnd,1)
  prevailingMean = matrix(0,TT-estimationEnd,1)
  
  predictors.model.lasso = matrix(0,TT-estimationEnd,10)


# Construct recursive forecasts
# You'll need to add some code that figures out which variables are
# selected using your preferred model selection method each time period




for (tt in seq(estimationEnd,TT-1)) {

  # say when tt <- 516 (1969-m12)
  print(TT-tt)    #display how many iterations are left to track progress, estimation should take several minutes
  # X[1:(tt-1),] is 1927-m01 until 1969-m11 all predictor variable values
  # y[2:tt] is 1927-m02 until 1969-m12 dependent variable values. Note a lag is created for the dependent variable
  # X[tt,] is 1969-m12 all predictor variable values
  # 2^10 =1024 different models are possible with 10 predictor variables
  # forecasts creates a matrix of 1024*1
  # criteria creates a matrix of 1024*2
  forecasts = Forecasts2K(X[1:(tt-1),],y[2:tt],X[tt,])$AllF     #constructs forecasts from every combination of predictors
  criteria = Forecasts2K(X[1:(tt-1),],y[2:tt],X[tt,])$Criteria  # holds AIC/Akaike and BIC values for each model 
  # Remember lower BIC values the better model
  minIndices = apply(criteria,2,which.min)                      # Apply BIC criteria and select the minimum from the different models created - for each iteration different model implies #of models based on # of predictor variables. 10 in this case
  cv = cv.glmnet(X[1:(tt-1),],y[2:tt],nfolds = 10)           #estimate a LASSO regression and keep one with lowest MSE using 10-fold cross-validation
  # coef(cv,s='lambda.min') will return a 11*1 matrix along with the coefficient values for each variables(10) + intercept corresponding to the lowest lambda value
  # !=0 check will return true /false if or not the coefficents are zero
  # -1 to get rid of the intercept
  IndexMinMSE = (as.vector(coef(cv,s='lambda.min'))!=0)[-1]
  bLasso = numeric(ncol(X)+1)
  # Training set will take only those predictor values which are True in IndexMinMSE(i.e corresponding to the lowest lambda value )
  trainset = cbind(rep(1,tt-1),X[1:(tt-1),IndexMinMSE])
  tempLasso = OLS(trainset,y[2:tt])$Beta
  
  # get bLasso coefficient
  bLasso[which(IndexMinMSE)+1] = tempLasso[-1]
  bLasso[1] = tempLasso[1]
  # get the predictors selected in each iteration
  predictors.model.lasso[tt-estimationEnd+1,] = IndexMinMSE
  
  # stepwise regression
  statsBack=step(lm(y[2:tt]~.,as.data.frame(X[1:(tt-1),])), direction="backward",test='F',trace = 0) #backward stepwise model selection
  statsForward=step(lm(y[2:tt]~1,data=as.data.frame(X[1:(tt-1),])), direction="forward",test='F',formula(lm(y[2:tt]~.,as.data.frame(X[1:(tt-1),]))),trace=0) #forward stepwise model selection
  IndexForward=colnames(X)%in%names(coef(statsForward))
  IndexBackward=colnames(X)%in%names(coef(statsBack))
  
  # ?lowest isn't it
  AICforecast[tt-estimationEnd+1,] = forecasts[minIndices[1]] #highest AIC forecast
  BICforecast[tt-estimationEnd+1,] = forecasts[minIndices[2]] #highest BIC forecast
  
  # 516th row predictor variables values
  now = data.frame(X)[tt,]
  forwardStepForecast[tt-estimationEnd+1] = predict(statsForward,newdata=now) #forward stepwise forecast
  backStepForecast[tt-estimationEnd+1] = predict(statsBack,newdata=now)       #backward stepwise forecast
  combinationForecast[tt-estimationEnd+1] = mean(forecasts)                   #equal-weighted average of all model forecasts
  lassoForecast[tt-estimationEnd+1] = sum(c(1,X[tt,])*bLasso)                 #min MSE Lasso forecast
  kitchenSink[tt-estimationEnd+1] = forecasts[length(forecasts)]              #kitchen sink model forecast
  prevailingMean[tt-estimationEnd+1] = mean(y[2:tt])                          #prevailing mean forecast
}

#Aggregate forecasts into matrix and construct forecast errors
forecasts = cbind(AICforecast, BICforecast, forwardStepForecast, backStepForecast, combinationForecast, lassoForecast, kitchenSink, prevailingMean)
forecastsErrors = y[seq(estimationEnd,TT-1)]-forecasts

#Plot model forecasts
timeVector = seq(1969+11/12,2015+11/12,length.out =  TT-estimationEnd)  #you need to change this
matplot(timeVector,forecasts,type = 'l',col=palette(),lty=1,xlab = 'Time',ylab = 'Percentage Points')
legend(min(timeVector),-0.01,modelLabels[-1],col=palette(),lty = rep(1,8))

#Compute model RMSEs (root mean squared errors)
RMSE = sqrt(colMeans(forecastsErrors^2))
RMSE
### Q1 : Looking at the Root mean Squared error, we deduce that Combination Forecast (which is the )equal-weighted average of all model forecasts) produces best forecast with least RMSE.
## However if we consider models by AIC or BIC,LASSO and stepwise selection methods , we can see LASSO having the edge with least RMSE among them.

## Here is the plot of actuals vs forecast for each of the methods .

for(i in (1:8)){
  #Plot model forecasts 
  timeVector = seq(1969+11/12,2015+11/12,length.out =  TT-estimationEnd)  #you need to change this
  forecast.actual <- cbind(y[seq(estimationEnd,TT-1)],forecasts[,i])
  matplot(timeVector,forecast.actual,type = 'l',col=palette(),lty=1,xlab = 'Time',ylab = 'Percentage Points')
  forecast.actual.legend <- c("actual",modelLabels[-1][i])
  legend(min(timeVector),-0.01,forecast.actual.legend,col=palette(),lty = rep(1,2))
}


### Q2: How often does your preferred model selection approach include different predictor variables?

predictor.variables <- cbind(colnames(datas[seq(1,nrow(datas)),c(6:15)]),round(apply(predictors.model.lasso,2,mean)*100,2))
predictor.variables
## book-to-market ratio is included in the model 30.43% of the time followed by net equity issues which is included 17.2% of the time
## default spread,dividend-price ratio and stock variance are less likely to play an important role in the forecasting of excess stock return over time as they have been included in only 0.18%,1.45% and 1.45% of the models 


### Q3: 
# LASSO model model produces better forecast with RMSE of 0.04431 when compared to the Prevailing mean model with RMSE of 0.04432

### Q4: 
# LASSO model model produces better forecast with RMSE of 0.04431 when compared to the Kitchen Sink model with RMSE of 0.04498


### Q5:
# Calculate economic performance of models (relative to risk free rate)
riskFree = datas[seq(1,nrow(datas)),17] #you need to change this
stockReturn = datas[seq(1,nrow(datas)),ncol(datas)-1] # you need to change this
economicReturns = matrix(0,nrow(forecasts),ncol(forecasts))

## Changed seq(1:estimationEnd) to seq(estimationEnd,TT-1)
for (ii in seq(1,ncol(forecasts))) {
  negativeReturns = forecasts[,ii]<=0
  economicReturns[,ii] = riskFree[seq(estimationEnd,TT-1)]*negativeReturns+stockReturn[seq(estimationEnd,TT-1)]*(!negativeReturns)
}

# Compute mean returns and Sharpe Ratios

meanReturns = colMeans(economicReturns)
meanReturns
sharpeRatios = (meanReturns-mean(riskFree[-seq(1,estimationEnd)]))/apply(economicReturns, 2, sd)
sharpeRatios

## Combination Forecast model generates highest mean Returns and Sharp Ratios.
## Among the three models discussed above Prevailing Mean model generates highest mean Returns and Sharp Ratios

## mean Returns of the three models discussed above
# prevailingMean	0.009225928
# lassoForecast,	0.008953489
# kitchenSink,	0.007790466
# 
# Sharp Ratios of the three models discussed above
# prevailingMean	0.11564359
# lassoForecast,	0.11286074
# kitchenSink,	  0.10419858

save.image(file = "assignment2.rda")
#save.image(file = "assignment2.rda") # use this command to save your results

### Q6: 6.	Discuss your findings - what are your conclusions regarding how easy it is to predict stock returns?