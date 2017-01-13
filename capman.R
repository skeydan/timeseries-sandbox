library(readr)
library(dplyr)
library(forecast)
library(ggplot2)
library(lubridate)
library(dlm)

dbcpu <- read_csv2('dbcpu.csv', col_names = c('rownum','tstamp','value'),
                  col_types= 'icn_', skip=1)
dbcpu <- dbcpu %>% mutate(tstamp = parse_datetime(tstamp))
dbcpu <- dbcpu %>% filter(tstamp > '2016-01-01')


############################ interventions ###################################

# 16 April, 2016 Hardware-Migration ivorapo02
dbcpu <- dbcpu %>% mutate(intervention1 = ifelse(tstamp > '2016-04-16', 1, 0))

# 6 June, 2016 Hardware-Migration ivorapo01
dbcpu <- dbcpu %>% mutate(intervention2 = ifelse(tstamp > '2016-06-06', 1, 0))

# 24 July, 2016 cpupower frequency-set --governor performance
dbcpu <- dbcpu %>% mutate(intervention3 = ifelse(tstamp > '2016-07-24', 1, 0))

dbtime <- read_csv2('dbtime.csv', col_names = c('rownum','tstamp','value'),
                   col_types= 'icn_', skip=1)
dbtime <- dbtime %>% mutate(tstamp = parse_datetime(tstamp))
head(dbtime)

iops <- read_csv2('iops.csv', col_names = c('rownum','tstamp','value'),
                   col_types= 'icn_', skip=1)
iops <- iops %>% mutate(tstamp = parse_datetime(tstamp))
head(iops)

size <- read_csv2('size.csv', col_names = c('rownum','tstamp','value'),
                   col_types= 'icn_', skip=1)
size <- size %>% mutate(tstamp = parse_datetime(tstamp))
head(size)


############################ plot ###################################

# "lm", "glm", "gam", "loess", "rlm"
ggplot(dbcpu, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbcpu") + geom_smooth(method = 'auto')
ggplot(dbcpu, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbcpu") + geom_smooth(method = 'loess', span = 0.1)

ggplot(dbtime, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbcpu") + geom_smooth(method = 'auto')
ggplot(dbtime, aes(tstamp, value)) + geom_point() + xlab("") + ylab("dbtime") + stat_smooth(method = 'loess', span = 0.1)

ggplot(iops, aes(tstamp, value)) + geom_point() + xlab("") + ylab("iops") + stat_smooth()
ggplot(size, aes(tstamp, value)) + geom_point() + xlab("") + ylab("size") + stat_smooth()

dbcpu_ts <- ts(dbcpu$value)

acf(dbcpu_ts)
pacf(dbcpu_ts)

plot(dbcpu_ts)
plot(diff(dbcpu_ts))
ndiffs(dbcpu_ts)

############################ ets ###################################

fit <- ets(dbcpu_ts)
summary(fit)
plot(fit)
plot(forecast(fit, h=90))

############################ auto.arima #############################

fit <- auto.arima(dbcpu_ts)
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
buildLocalLevel <- function(psi) dlmModPoly(order = 1, dV = psi[1], dW = psi[2])
mleOut <- dlmMLE(dbcpu_ts, parm = c(0.2, 120), build = buildLocalLevel, lower = c(1e-7, 0))
# Checking that the MLE estimates has converged
mleOut$convergence
# Constructing the fitted model
LocalLevelmod <- buildLocalLevel(psi = mleOut$par)
# observation variance
drop(V(LocalLevelmod))
# system variance
drop(W(LocalLevelmod))

