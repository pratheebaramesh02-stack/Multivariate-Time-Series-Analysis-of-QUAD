#TIME SERIES ANALYSIS OF THE STOCK INDICES OF THE QUAD COUNTRIES

#Our objective is to find out any long run relationship among the stock indices of the
#4 counties - India, USA, Japan and Australia. The idea is to look for cointegrating relationships
#for short run and long run dynamics in the model and find out causalities and impacts of shocks on
#the movement of these indices

#we start by clearing variables in workspace
rm(list = ls())

#installing packages
install.packages("ggplot2")
install.packages("zoo")
install.packages("readxl")
install.packages("tseries")
install.packages("fpp2")
install.packages('urca')
install.packages("vars")
install.packages("tsDyn")
install.packages("lmtest")

#importing packages
library(ggplot2)
library(zoo)
library(readxl)
library(tseries)
library(fpp2)
library(urca)
library(vars)
library(tsDyn)
library(lmtest)


#importing the excel sheet
multi <- read_excel("C:/Users/91994/Desktop/amfe_project/STIND.xlsx")
View(multi)

#making all the series into time series objects and getting their log
lnsei <- log(ts(multi[,2],start = c(2013,1),frequency = 52))
lgspc <- log(ts(multi[,3],start = c(2013,1),frequency = 52))
ln225 <- log(ts(multi[,4],start = c(2013,1),frequency = 52))
laxjo <- log(ts(multi[,5],start = c(2013,1),frequency = 52))

#making a multi-variate time series with all 4
data <- ts.union(lnsei,lgspc,ln225,laxjo)

#Plotting the series to check for stationarity
plot(data)

#Plotting them together for a better comparison
autoplot(data) +
  ggtitle("Time Series Plot of the indices - before differencing")

#Formally checking for a unit root with adf test
adf.test(lnsei)
adf.test(lgspc)
adf.test(ln225)
adf.test(laxjo)

#Formally checking for a unit root with kpss test
kpss.test(lnsei)
kpss.test(lgspc)
kpss.test(ln225)
kpss.test(laxjo)

#The results of adf and kpss suggest that the data is not stationary at level

#First-differencing the data
dlnsei <- diff(lnsei)
dlgspc <- diff(lgspc)
dln225 <- diff(ln225)
dlaxjo <- diff(laxjo)

#making a multi-variate time series with all 4
data1 <- ts.union(dlnsei,dlgspc,dln225,dlaxjo)

#Plotting the series to check for stationarity
plot(data1)

#Plotting them together for a better comparison
autoplot(data1) +
  ggtitle("Time Series Plot of the indices - after differencing")


#Formally checking for stationarity with adf test
adf.test(dlnsei)
adf.test(dlgspc)
adf.test(dln225)
adf.test(dlaxjo)

#Formally checking for stationarity with kpss test
kpss.test(dlnsei)
kpss.test(dlgspc)
kpss.test(dln225)
kpss.test(dlaxjo)

#The results of adf and kpss suggest that the data has become stationary at
#first-difference, indicating we are working with I(1) variables

#Choosing lag length for the JJ test
VARselect(data1, lag.max = 10, type = "none")$selection

#AIC and FPE shows 2, whereas BIC and HQIC shows 1. Since the JJ test in R 
#requires a minimum lag length of 2, we go ahead with 2

#Johansen cointegration test with trace statistic
result <- ca.jo(data, type = "trace", K=2, ecdet = "trend",spec="transitory")
summary(result)

#JJ test results show 2 cointegrating relationships. The next step is a VECM

#VECM (Interpretation of the same in the report)
vecm <- VECM(data.frame(lnsei,ln225,laxjo,lgspc), lag = 2, r = 2, include = "both", estim = "ML")
summary(vecm)

#We convert the vecm to a var to conduct diagnostic checks and other analysis.

#VECM to VAR
var <- vec2var(result, r = 2)

#Checking for autocorrelaion in residuals
serial.test(var, lags.pt = 18)

#We fail to reject H0 of no autocorrelation at 1% and 5% levels

#Granger causality for the variables
grangertest(dlgspc ~ dlnsei)
grangertest(dlgspc ~ dln225)
grangertest(dlgspc ~ dlaxjo)
grangertest(dlnsei ~ dlaxjo)
grangertest(dlnsei ~ dlgspc) #significant - GSPC granger causes NSEI
grangertest(dlnsei ~ dln225)
grangertest(dlaxjo ~ dln225) 
grangertest(dlaxjo ~ dlnsei)
grangertest(dlaxjo ~ dlgspc)
grangertest(dln225 ~ dlaxjo) 
grangertest(dln225 ~ dlnsei)
grangertest(dln225 ~ dlgspc) #significant  - GSPC granger causes N225

#Impulse Response Functions for the model

#Plotting the Implulse of GSPC on other indices
ir1 <- irf(var, n.ahead = 100, impulse = "lgspc", 
          response = c("lnsei", "ln225", "laxjo"),
          ortho = FALSE, runs = 500)
plot(ir1)

#Plotting the Implulse of NSEI on other indices
ir2 <- irf(var, n.ahead = 100, impulse = "lnsei",
           response = c("lgspc", "ln225", "laxjo"),
           ortho = FALSE, runs = 500)
plot(ir2)

#Plotting the Implulse of N225 on other indices
ir3 <- irf(var, n.ahead = 100, impulse = "ln225",
           response = c("lnsei", "lgspc", "laxjo"),
           ortho = FALSE, runs = 500)
plot(ir3)

#Plotting the Implulse of AXJO on other indices
ir4 <- irf(var, n.ahead = 100, impulse = "laxjo",
           response = c("lnsei", "ln225", "lgspc"),
           ortho = FALSE, runs = 500)
plot(ir4)

#Error Variance Decomposition for the model
fevd_result <- fevd(var, n.ahead = 10)

#Printing out the FEVD numerically
print(fevd_result)

#plotting the bars to check visually
plot(fevd_result)


###############################################################################
