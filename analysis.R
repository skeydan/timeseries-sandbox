library(dplyr)
library(readr)
library(ggplot2)
library(forecast)
library(stlplus)


id <- substr(filename,1,3)
reason <- 7
tank_nbr <- 3
measurement_error_threshold <- -10
test_set_num_days <- 90

df <- read_csv(filename, col_types = cols(`Mess-Dat` = col_date('%d/%m/%y')))

df <- df %>% filter(MessGrund == reason & TankNr == tank_nbr)
df <- df %>% select(`Mess-Dat`, `Mess-Zeit`, InhaltLiter)
df <- df %>% rename(day = `Mess-Dat`, time = `Mess-Zeit`, filling_height = InhaltLiter)
tail(df)
nrow(df)



# there are several measurements (not just at the beginning!) where obviously there was a problem and filling_height is 0.
# if we don't set them to NA, they will show up as crazy outliers in the diffs calculation later.
df <- df %>% mutate(filling_height = ifelse(filling_height == 0, NA, filling_height))
head(df)
df %>% filter(is.na(filling_height)) %>% nrow()

# filling_height minus previous filling_height is the amount that's gone between 2 measuring points
# diff normally should be >= 0
df <- df %>% mutate(prev_height = lag(filling_height)) %>% mutate(diff = prev_height - filling_height) 
tail(df)
df %>% filter(is.na(diff)) %>% nrow()

# small negative diffs probably are measurement errors, we set them to 0
# big negative diffs are refills, we mark them as such but set them to 0 for calculation
# this introduces an error (we miss how much fuel was taken up in the same time interval)
df <- df %>% mutate(diff_cleaned = ifelse(diff < 0, 0, diff),
                    refill = ifelse(diff < measurement_error_threshold, 1, 0)
)
tail(df)

# sanity check
# we should see that in this selection
## diff_cleaned is 0 everywhere
## refill is 1 only when diff is smaller than -10
df %>% filter(diff < -5) %>% tail(20)


# NA handling step 1: remove all missing values at the start of the measurement period
nrow(df)
first_not_na <- df %>% filter(!is.na(diff_cleaned)) %>% head(1) %>% select(day,time)
df <- filter(df, day > first_not_na$day)
nrow(df)

# NA handling step 2: 
# how many missing value per day?
df %>% filter(is.na(diff_cleaned)) %>% group_by(day) %>% summarise(count = n())
# example output
#1 2015-04-23     2
#2 2015-06-15    11
#3 2015-06-16    24
#4 2015-06-17    24
#5 2015-06-18    16
#6 2016-02-18     7
#7 2016-03-16    12
#8 2016-03-17    11
#9 2016-10-04     9

# thus for simplicity, we always insert mean value instead of interpolating
mean_diff <-(df %>% summarize(mean_diff = mean(diff_cleaned, na.rm = TRUE)))$mean_diff
mean_diff
df <- df %>% mutate(diff_cleaned = ifelse(is.na(diff_cleaned), mean_diff, diff_cleaned))
# this must be 0 now!
df %>% filter(is.na(diff_cleaned)) %>% group_by(day) %>% summarise(count = n())


# aggregate by day
df_day <- df %>%  group_by(day) %>% summarize(agg_diff = sum(diff_cleaned), refill = sum(refill, na.rm = TRUE))
df_day %>% filter(is.na(agg_diff))
head(df_day)
tail(df_day)
nrow(df_day)




#################################################################
###                             plot                          ### 
#################################################################

# there should not be any extreme outliers here now, if there are check why!
ggplot(df_day, aes(day,agg_diff)) + geom_line()
# zoom in
start_date <- as.Date('2016-01-01')
end_date <- as.Date('2017-01-01')
ggplot(df_day, aes(day,agg_diff)) + geom_line() + coord_cartesian(xlim = c(start_date, end_date))




#################################################################
###              split into test and training sets            ### 
#################################################################


# split into test and training sets
df_day_train <- head(df_day, -test_set_num_days)
df_day_test <- tail(df_day, test_set_num_days)

write_csv(df_day_train, paste0(paste(id, reason, tank_nbr, sep = '_'),'.train.csv'))
write_csv(df_day_test, paste0(paste(id, reason, tank_nbr, sep = '_'),'.test.csv'))

#################################################################
###                       arima                               ### 
#################################################################

arima_fit <- auto.arima(df_day_train$agg_diff, trace = TRUE)
summary(arima_fit)

fc <- forecast(arima_fit, h = nrow(df_day_test))
summary(fc)

# fc$fitted gives the one-step-ahead forecasts
# fc$residuals is the one-step-ahead forecast error
fc$x - fc$fitted - fc$residuals

plot(fc)


ets_fit <- ets(df_day_train$agg_diff)
summary(ets_fit)
plot(forecast(ets_fit, h= nrow(df_day_test)))

#################################################################
###            use ts objects to specify seasonality          ### 
#################################################################

ts_day_train <- ts(df_day_train$agg_diff, deltat = 1/7)

plot(stlplus(ts_day_train, s.window = 'periodic'))

arima_fit_ts <- auto.arima(ts_day_train, trace = TRUE)
summary(arima_fit_ts)
fc_ts <- forecast(arima_fit_ts, h = nrow(df_day_test))
fc_ts <- forecast(arima_fit_ts, h = 3)

summary(fc_ts)

plot(fc_ts)

ets_fit_ts <- ets(ts_day_train)
summary(ets_fit_ts)
plot(forecast(ets_fit_ts, h= nrow(df_day_test)))

