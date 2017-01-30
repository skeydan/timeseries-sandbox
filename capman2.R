library(readr)
library(dplyr)
library(forecast)
library(ggplot2)
library(lubridate)
library(dlm)

# forecast ab 14.1.

dbcpu1 <- read_csv2('dbcpu_daily_with_forecast.csv', col_names = c('rownum','tstamp','value'),
                  col_types= 'icn_', skip=1)
dbcpu1 <- dbcpu1 %>% mutate(tstamp = parse_datetime(tstamp))
dbcpu1 <- dbcpu1 %>% mutate(actual = if_else(tstamp < '2017-01-14', value, as.double(NA)),
                            pred = if_else(tstamp < '2017-01-14', as.double(NA), value))
head(dbcpu1)
tail(dbcpu1)

dbcpu2 <- read_csv2('dbcpu_hourly_with_prediction.csv', col_names = c('rownum','tstamp','value'),
                    col_types= 'icn_', skip=1)
dbcpu2 <- dbcpu2 %>% mutate(tstamp = parse_datetime(tstamp))
dbcpu2 <- dbcpu2 %>% mutate(actual = if_else(tstamp < '2017-01-14', value, as.double(NA)),
                            pred = if_else(tstamp < '2017-01-14', as.double(NA), value))
head(dbcpu2)
tail(dbcpu2)

dbcpu <- dbcpu2
#dbcpu_ts <- ts(dbcpu$actual)
#dbcpu_ts <- ts(dbcpu$actual, frequency=365)
dbcpu_ts <- ts(dbcpu$actual, frequency=8)
############################ plot ###################################

# "lm", "glm", "gam", "loess", "rlm"
ggplot(dbcpu, aes(x = tstamp)) + geom_point(aes(y = actual), color='blue') + geom_point(aes(y = pred), color = 'green') + 
  geom_smooth(aes(y = actual))
  #geom_smooth(aes(y = actual), method = 'loess', span = 0.01)


#acf(dbcpu_ts)
#pacf(dbcpu_ts)
#plot(dbcpu_ts)
#plot(diff(dbcpu_ts))
#ndiffs(dbcpu_ts)

############################ ets ###################################

#fit <- ets(dbcpu_ts)
#summary(fit)
#plot(fit)
#plot(forecast(fit, h=90))

############################ auto.arima #############################

#default
fit <- auto.arima(dbcpu_ts)
summary(fit)
#fit <- auto.arima(dbcpu_ts, max.order = 10, stepwise=FALSE)
plot(forecast(fit, h=90))


# fourier
K = 4
h = 90
fit <- auto.arima(dbcpu_ts, seasonal=FALSE, xreg=fourier(dbcpu_ts, K=K))
summary(fit)
plot(forecast(fit, h=h, xreg = fourier(dbcpu_ts, K=K, h=h)))


############################ kalman ###################################

ts <- dbcpu_ts
#ts <- iops_ts
#ts <- size_ts

#  Random walk plus noise model (polynomial model of order one)

# dlmModPoly(order = 2, dV = 1, dW = c(rep(0, order - 1), 1), m0 = rep(0, order), C0 = 1e+07 * diag(nrow = order)) 
# dV = variance of the observation noise
# dW = diagonal elements of the variance matrix of the system noise
# build a local level model
build_model <- function(params) dlmModPoly(order = 1, dV = params[1], dW = params[2])

# dlmMLE(y, parm, build, method = "L-BFGS-B", ..., debug = FALSE) 
# parm: vector of initial values for optimization
# build: function from vector of same length as parm to a dlm object
estim <- dlmMLE(ts, parm = c(0,0), build = build_model)

# Checking that the MLE estimates has converged
estim$convergence

# Constructing the fitted model
fitted <- build_model(params = estim$par)
# observation variance
fitted$V
# system variance
fitted$W

filtered<- dlmFilter(y = ts, mod = fitted)
smoothed <- dlmSmooth(y = ts, mod = fitted)
# If you do end up with a time-variant model, you'll want to fill out your input data with NA's and let 
# the dlmFilter fill in the NA's for you (a poor man's forecast), since dlmForecast does not work with time-varying
# parameters.
forecast_obj <- dlmForecast(fitted, nAhead = 30)
forecast_obj$f

head(filtered$m)
head(smoothed$s)
head(ts)
tail(filtered$f)
tail(smoothed$s)
tail(ts)

comp_df <- data_frame(x=index(ts), orig=ts, f = filtered$m[-1], s = smoothed$s[-1])
ggplot(comp_df, aes(x)) + 
  geom_point(aes(y=orig), colour='black', size=0.2) +
  geom_line(aes(y=f), colour='red', size=0.2) + 
  geom_line(aes(y=s), colour='green', size=0.2)

