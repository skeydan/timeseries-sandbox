library(readr)
library(dplyr)
library(forecast)
library(ggplot2)
library(lubridate)
library(dlm)

dbcpu <- read_csv2('dbcpu.csv', col_names = c('rownum','tstamp','value'),
                  col_types= 'icn_', skip=1)
dbcpu <- dbcpu %>% mutate(tstamp = parse_datetime(tstamp))
#dbcpu <- dbcpu %>% filter(tstamp > '2016-01-01')


############################ interventions ###################################

# 16 April, 2016 Hardware-Migration ivorapo02
dbcpu <- dbcpu %>% mutate(intervention1 = ifelse(tstamp > '2016-04-16', 1, 0))

# 6 June, 2016 Hardware-Migration ivorapo01
dbcpu <- dbcpu %>% mutate(intervention2 = ifelse(tstamp > '2016-06-06', 1, 0))

# 24 July, 2016 cpupower frequency-set --governor performance
dbcpu <- dbcpu %>% mutate(intervention3 = ifelse(tstamp > '2016-07-24', 1, 0))

head(dbcpu)

dbtime <- read_csv2('dbtime.csv', col_names = c('rownum','tstamp','value'),
                   col_types= 'icn_', skip=1)
dbtime <- dbtime %>% mutate(tstamp = parse_datetime(tstamp))
#dbtime <- dbtime %>% filter(tstamp > '2016-01-01')
head(dbtime)

iops <- read_csv2('iops.csv', col_names = c('rownum','tstamp','value'),
                   col_types= 'icn_', skip=1)
iops <- iops %>% mutate(tstamp = parse_datetime(tstamp))
#iops <- iops %>% filter(tstamp > '2016-01-01')
head(iops)

size <- read_csv2('size.csv', col_names = c('rownum','tstamp','value'),
                   col_types= 'icn_', skip=1)
size <- size %>% mutate(tstamp = parse_datetime(tstamp))
#size <- size %>% filter(tstamp > '2016-01-01')
head(size)


############################ plot ###################################

# "lm", "glm", "gam", "loess", "rlm"
ggplot(dbcpu, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbcpu") + geom_smooth(method = 'auto')
ggplot(dbcpu, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbcpu") + geom_smooth(method = 'loess', span = 0.01)

ggplot(dbtime, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbcpu") + geom_smooth(method = 'auto')
ggplot(dbtime, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbtime") + stat_smooth(method = 'loess', span = 0.1)

ggplot(iops, aes(tstamp, value)) + geom_point() + xlab("") + ylab("iops") + stat_smooth()
ggplot(size, aes(tstamp, value)) + geom_point() + xlab("") + ylab("size") + stat_smooth()

dbcpu_ts <- ts(dbcpu$value)
dbtime_ts <- ts(dbtime$value)
iops_ts <- ts(iops$value)
size_ts <- ts(size$value)

acf(dbcpu_ts)
acf(dbtime_ts)
acf(iops_ts)
acf(size_ts)

pacf(dbcpu_ts)
pacf(dbtime_ts)
pacf(iops_ts)
pacf(size_ts)

plot(dbcpu_ts)
plot(diff(dbcpu_ts))

ndiffs(dbcpu_ts)
ndiffs(dbtime_ts)
ndiffs(iops_ts)
ndiffs(size_ts)


############################ ets ###################################

fit <- ets(dbcpu_ts)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

fit <- ets(dbtime_ts)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

fit <- ets(iops_ts)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

fit <- ets(size_ts)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

############################ auto.arima #############################

fit <- auto.arima(dbcpu_ts, max.order = 20, stepwise=FALSE)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

fit <- auto.arima(dbtime_ts, max.order = 20, stepwise=FALSE)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

fit <- auto.arima(iops_ts, max.order = 20, stepwise=FALSE)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

fit <- auto.arima(size_ts, max.order = 20, stepwise=FALSE)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))


############################ regression with arima errors #################

reg <- cbind(dbcpu$intervention1, dbcpu$intervention2, dbcpu$intervention3)
fit <- auto.arima(dbcpu_ts, xreg = reg)
summary(fit)
plot(fit)
plot(forecast(fit, xreg=cbind(rep(1,90), rep(1,90), rep(1,90)), h=90))


############################ kalman ###################################

#  Random walk plus noise model (polynomial model of order one)

# dlmModPoly(order = 2, dV = 1, dW = c(rep(0, order - 1), 1), m0 = rep(0, order), C0 = 1e+07 * diag(nrow = order)) 
# dV = variance of the observation noise
# dW = diagonal elements of the variance matrix of the system noise
# build a local level model
build_model <- function(params) dlmModPoly(order = 1, dV = params[1], dW = params[2])

# dlmMLE(y, parm, build, method = "L-BFGS-B", ..., debug = FALSE) 
# parm: vector of initial values for optimization
# build: function from vector of same length as parm to a dlm object
estim <- dlmMLE(dbcpu_ts, parm = c(0,0), build = build_model)

# Checking that the MLE estimates has converged
estim$convergence

# Constructing the fitted model
fitted <- build_model(params = estim$par)
# observation variance
fitted$V
# system variance
fitted$W

filtered<- dlmFilter(y = dbcpu_ts, mod = fitted)
smoothed <- dlmSmooth(y = dbcpu_ts, mod = fitted)
# If you do end up with a time-variant model, you'll want to fill out your input data with NA's and let 
# the dlmFilter fill in the NA's for you (a poor man's forecast), since dlmForecast does not work with time-varying
# parameters.
forecast_obj <- dlmForecast(fitted, nAhead = 30)
forecast_obj$f

head(filtered$m)
head(smoothed$s)
head(dbcpu_ts)
tail(filtered$f)
tail(smoothed$s)
tail(dbcpu_ts)

comp_df <- data_frame(x=index(dbcpu_ts), orig=dbcpu_ts, f = filtered$m[-1], s = smoothed$s[-1])
ggplot(comp_df, aes(x)) + 
  geom_point(aes(y=orig), colour='black', size=0.2) +
  geom_line(aes(y=f), colour='red', size=0.2) + 
  geom_line(aes(y=s), colour='green', size=0.2)

