library(dlm)
library(PerformanceAnalytics)

# from: # http://lalas.github.io/quantitativeThoughts/r/2014/09/01/dlmTutorial.html


###########################################################################
#                Time varying univariate linear regression                #
###########################################################################

data(managers)
# extract HAM1 and SP500 excess returns
HAM1 = 100*(managers[,"HAM1", drop=FALSE] - managers[,"US 3m TR", drop=FALSE])
sp500 = 100*(managers[,"SP500 TR", drop=FALSE] - managers[,"US 3m TR",drop=FALSE])
colnames(sp500) = "SP500"

# Specifying a set model parameters
# example
#s2_obs = 1      # Variance of observations
#s2_alpha = 0.01 # Variance of the alpha regression parameter
#s2_beta = 0.01  # Variance of the beta regression parameter

# Construct a regression model
#tvp.dlm = dlmModReg(X=sp500, addInt=TRUE, dV=s2_obs, dW=c(s2_alpha, s2_beta))

# looking at the various components
# FF = observation matrix, GG = transition matrix
# V = covariance matrix of observation errors, W = cov matrix of transition errors
# m0 = initial mean of state, C0 = initial covariance of state
#tvp.dlm[c("FF","V","GG","W","m0","C0")]

# 
# X: used to store all the time-varying elements of the model
# JFF, JV, JGG, JW:  indicator matrices whose entries signal whether an element of the corresponding model matrix is time-varying and, in case it is, where to retrieve its values in the matrix X.
#tvp.dlm[c("JFF","JV","JGG","JW")]
#str(tvp.dlm["X"])

# for comparison: OLS
#ols.fit = lm(HAM1 ~ sp500)
#summary(ols.fit)


# s_obs, s_alpa and s_beta estimation
start.vals = c(0,0,0)
names(start.vals) = c("lns2_obs", "lns2_alpha", "lns2_beta")

# function to build Time Varying Parameter state space model
buildTVP <- function(parm, x.mat){
  parm <- exp(parm)
  return( dlmModReg(X=x.mat, dV=parm[1], dW=c(parm[2], parm[3])) )
}

# Estimate the model
TVP.mle = dlmMLE(y=HAM1, parm=start.vals, x.mat=sp500, build=buildTVP, hessian=T)

# get sd estimates
se2 <- sqrt(exp(TVP.mle$par))
names(se2) = c("s_obs", "s_alpha", "s_beta")
sqrt(se2)

# build model with estimated parameters
TVP.dlm <- buildTVP(TVP.mle$par, sp500)
TVP.dlm

# filtering
TVP.f <- dlmFilter(y = HAM1, mod = TVP.dlm)
class(TVP.f)
names(TVP.f)

# smoothing
TVP.s <- dlmSmooth(TVP.f)
class(TVP.s)

alpha.s = xts(TVP.s$s[-1,1,drop=FALSE], as.Date(rownames(TVP.s$s[-1,])))
beta.s  = xts(TVP.s$s[-1,2,drop=FALSE], as.Date(rownames(TVP.s$s[-1,])))
colnames(alpha.s) = "alpha"
colnames(beta.s)  = "beta"

# extract std errors - dlmSvd2var gives list of MSE matrices
mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)
se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
se.xts = xts(se.mat[-1, ], index(beta.s))
colnames(se.xts) = c("alpha", "beta")
a.u = alpha.s + 1.96*se.xts[, "alpha"]
a.l = alpha.s - 1.96*se.xts[, "alpha"]
b.u = beta.s  + 1.96*se.xts[, "beta"]
b.l = beta.s  - 1.96*se.xts[, "beta"]

# plot smoothed estimates with +/- 2*SE bands
chart.TimeSeries(cbind(alpha.s, a.l, a.u), main="Smoothed estimates of alpha", ylim=c(0,1),
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(alpha),xlab="")


chart.TimeSeries(cbind(beta.s, b.l, b.u), main="Smoothed estimates of beta",
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(beta),xlab="")

library(ggplot2, warn.conflicts = FALSE)
alpha.df <- data.frame(dateTime = index(se.xts), alpha = alpha.s, upr = a.u, lwr = a.l)
names(alpha.df) <- c("dateTime", "alpha", "upr", "lwr")
beta.df  <- data.frame(dateTime = index(se.xts), beta = beta.s, upr = b.u, lwr = b.l)
names(beta.df) <- c("dateTime", "beta", "upr", "lwr")

## Plotting alpha
ggplot(data = alpha.df, aes(dateTime, alpha) ) + geom_point () + geom_line() +
  geom_ribbon(data=alpha.df, aes(ymin=lwr,ymax=upr), alpha=0.3) + 
  labs(x = "year", y = expression(alpha), title = expression(paste("State Space Values of ", alpha, " over Time")))

## Plotting beta
ggplot(data = beta.df, aes(dateTime, beta) ) + geom_point (data = beta.df, aes(dateTime, beta) ) +
  geom_line() + geom_ribbon(data=beta.df , aes(ymin=lwr,ymax=upr), alpha=0.3) +
  labs(x = "year", y = expression(beta), title = expression(paste("State Space Values of ", beta, " over Time")))

# Construct add 10 missing values to end of sample
new.xts = xts(rep(NA, 10), seq.Date(from=end(HAM1), by="months", length.out=11)[-1])
# Add this NA data to the original y (HAM1) series
HAM1.ext = merge(HAM1, new.xts)[,1]
# Filter extended y (HAM1) series
TVP.ext.f = dlmFilter(HAM1.ext, TVP.dlm)
# extract h-step ahead forecasts of state vector
TVP.ext.f$m[as.character(index(new.xts)),]

TVP.res <- residuals(TVP.f, sd = FALSE)
# Q-Q plot
qqnorm(TVP.res)
qqline(TVP.res)

tsdiag(TVP.f)



###########################################################################
#                           Random walk with noise                        #
###########################################################################


data(Nile)
Nile.df <- data.frame(year = index(Nile), y = as.numeric(Nile))
qplot(y = y, x = year, data = Nile.df, geom = 'line', ylab = 'Nile water level', xlab = 'Year',
      main = "Measurements of the annual flow of \n the river Nile at Ashwan 1871-1970")

# Creating models -- the variance V is the same in mod1 and mod2; but 
# the signal variance is 10 time larger in mod2 than mod1
mod1 <- dlmModPoly(order = 1, dV = 15100, dW = 0.5 * 1468)
mod2 <- dlmModPoly(order = 1, dV = 15100, dW = 5 * 1468)
# Creating filter data
NileFilt_1 <- dlmFilter(Nile, mod1)
NileFilt_2 <- dlmFilter(Nile, mod2)
# Creating df to contain data to plot
Nile.df <- data.frame(year = index(Nile), Orig_data = as.numeric(Nile), 
                      Filt_mod_1 = as.numeric(NileFilt_1$m[-1]), Filt_mod_2 = as.numeric(NileFilt_2$m[-1]))

myColor <- c('green','blue','violet')
p <- ggplot(data = Nile.df) 
p <- p + geom_point(aes(x = year, y = Orig_data), size=0.5, colour= "black", shape = 21, fill = myColor[1])
p <- p + geom_line(aes(x = year, y = Orig_data,  colour = "Orig_data") , size = 0.5)
p <- p + geom_line(aes(x = year, y = Filt_mod_1, colour = "Filt_mod_1"), linetype="dotdash")
p <- p + geom_line(aes(x = year, y = Filt_mod_2, colour = "Filt_mod_2"), linetype="dotdash")
p <- p + labs(x = "year", y = "water Level", title = "Nile River water level for \n 2 different signal-to-noise ratios")
p <- p + scale_colour_manual("", breaks = c("Orig_data", "Filt_mod_1", "Filt_mod_2"),
                             labels= c("Org Data", "Filterd Data 1: SNR = x", "Filterd Data 2: SNR = 10x"),
                             values = myColor[c(2,3,1)])
p <- p + theme(legend.position="bottom")
print(p)


### Model 1: assuming constant variances V and W

# build a local level model
buildLocalLevel <- function(psi) dlmModPoly(order = 1, dV = psi[1], dW = psi[2])
mleOut <- dlmMLE(Nile, parm = c(0.2, 120), build = buildLocalLevel, lower = c(1e-7, 0))
# Checking that the MLE estimates has converged
mleOut$convergence
# Constructing the fitted model
LocalLevelmod <- buildLocalLevel(psi = mleOut$par)
# observation variance
drop(V(LocalLevelmod))
# system variance
drop(W(LocalLevelmod))


### Model 2: let W change over time (in 1899)

# Model Construction
buildDamEffect <- function(psi) {
  mod <- dlmModPoly(1, dV = psi[1], C0 = 1e8)
  # Creating the X matrix for the model -- For more info see Time Varying DLM section above
  X(mod) <- matrix(psi[2], nr = length(Nile))
  X(mod)[time(Nile) == 1899] <- psi[3]
  # Tell R that the values of W_t at any time are to be found in the first column of the matrix X
  JW(mod) <- 1
  return(mod)
}
# Model estimation through MLE
mleDamEffect <- dlmMLE(Nile, parm = c(0.2, 120, 20), build = buildDamEffect, lower = c(1e-7, 0, 0))
# Verify convergence
mleDamEffect$conv
# Construct the final model
damEffect <- buildDamEffect(psi = mleDamEffect$par)


### Model 3: linear model

# Model Construction
buildLinearTrend <- function(psi) dlmModPoly(2, dV = psi[1], dW = psi[2:3], C0 = diag(1e8, 2))
# Model Estimation
mleLinearTrend <- dlmMLE(Nile, parm = c(0.2, 120, 20),build = buildLinearTrend, lower = c(1e-7, 0, 0))
# Checking convergence
mleLinearTrend$conv
# Construct Final Model
linearTrend <- buildLinearTrend(psi = mleLinearTrend$par)

# MLE results checking
# We validate the results returned by the MLE method and gain confidence in the results by examining:
library(numDeriv)
# Local Level model
hs_localLevel <- hessian(function(x) dlmLL(Nile, buildLocalLevel(x)), mleOut$par)
all(eigen(hs_localLevel, only.values = TRUE)$values > 0) # positive definite?
# Damn Effect model
hs_damnEffect <- hessian(function(x) dlmLL(Nile, buildDamEffect(x)), mleDamEffect$par)
all(eigen(hs_damnEffect, only.values = TRUE)$values > 0) # positive definite?
# Linear Trend model
hs_linearTrend <- hessian(function(x) dlmLL(Nile, buildLinearTrend(x)), mleLinearTrend$par)
all(eigen(hs_linearTrend, only.values = TRUE)$values > 0) # positive definite?

### Models Comparison

#Model selection for DLMs is usually based on either of the following criteria:
  
# Forecasting accuracy – such as Mean square error (MSE), Mean absolute deviation (MAD), Mean absolute percentage error (MAPE).
# Information criteria – such as AIC, BIC
# Bayes factors and posterior model probabilities (in a Bayesian setting)
# If simulation is used, as it is the case in a Bayesian setting, then averages are calculated after discarding the burn-in samples (Petris, 2011).

# Creating variable to hold the results
MSE <- MAD <- MAPE <- U <- logLik <- N <- AIC <- c()
# Calculating the filtered series for each model
LocalLevelmod_filtered <- dlmFilter(Nile, LocalLevelmod)
damEffect_filtered     <- dlmFilter(Nile, damEffect)
linearTrend_filtered   <- dlmFilter(Nile, linearTrend)
# Calculating the residuals
LocalLevel_resid  <- residuals(LocalLevelmod_filtered, type = "raw", sd = FALSE)
damEffect_resid   <- residuals(damEffect_filtered, type = "raw", sd = FALSE)
linearTrend_resid <- residuals(linearTrend_filtered , type = "raw", sd = FALSE)
# If sampling was obtained through simulation then we would remove the burn-in samples as in the next line
# linearTrend_resid <- tail(linearTrend_resid, -burn_in)
#
# Calculating statistics for different models:
# 1 LocalLevelmod
MSE["Local Level"] <- mean(LocalLevel_resid ^2)
MAD["Local Level"] <- mean(abs(LocalLevel_resid ))
MAPE["Local Level"] <- mean(abs(LocalLevel_resid) / as.numeric(Nile))
logLik["Local Level"] <- -mleOut$value
N["Local Level"] <- length(mleOut$par)
# 2 Dam Effect
MSE["Damn Effect"] <- mean(damEffect_resid^2)
MAD["Damn Effect"] <- mean(abs(damEffect_resid))
MAPE["Damn Effect"] <- mean(abs(damEffect_resid) / as.numeric(Nile))
logLik["Damn Effect"] <- -mleDamEffect$value
N["Damn Effect"] <- length(mleDamEffect$par)
# 3 linear trend
MSE["linear trend"] <- mean(linearTrend_resid^2)
MAD["linear trend"] <- mean(abs(linearTrend_resid))
MAPE["linear trend"] <- mean(abs(linearTrend_resid) / as.numeric(Nile))
logLik["linear trend"] <- -mleLinearTrend$value
N["linear trend"] <- length(mleLinearTrend$par)
# Calculating AIC and BIC- vectorized, for all models at once
AIC <- -2 * (logLik - N)
BIC <- -2 * logLik + N * log(length(Nile))
# Building a dataframe to store the results
results <- data.frame(MSE = MSE, MAD = MAD, MAPE = MAPE, logLik = logLik, AIC = AIC, BIC = BIC, NumParameter = N)

# Producing a table
library(knitr)
kable(results, digits=2, align = 'c')



###########################################################################
#                Generate random model                                    #
###########################################################################

# Random model
# default is constant over time
ndim_obs <- 1
ndim_state <- 2
nobs <-10
# either all elements of FF are time-varying, or none are
JFF <- FALSE
JV <- FALSE
JGG <- FALSE
JW <- FALSE
r <- dlmRandom(ndim_obs, ndim_state, nobs)

# the model
r$mod
# simulated observations
r$theta
# observations
r$y

# forecast from model
# only for constant models
f <- dlmForecast(r$mod, nAhead = 3, sampleNew = 100)
# expected values of future states
f$a
# variances of future states
f$R
# expected values of future observations
f$f
# variances of future observations
f$Q
# sample states and observations from the forecast distribution
f$newStates
f$newObs
