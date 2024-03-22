## ----------------------Title--------------------------------------------------
#   Statistical analysis

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

# Script description: Analysis of whether growth parameters of packaged spinach from Arizona
# or California are significantly different 

#Null hypothesis: Growth parameters of baby spinach is not significantly different
#                 by state of cultivation

#Alternative hypothesis: Growth parameters of baby spinach is not significantly different
#                 by state of cultivation

## ----------------------Packages-----------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(reshape2); library(biogrowth); library(lmerTest)

## ----------------------Reading in Data----------------------------------------

param <- read.csv("output/primary_growth_model_parameters.csv", header = TRUE)
metadata <- read.csv("data/wrangled/wrangled_metadata.csv", header = TRUE)
harvest <- read.csv("data/wrangled/harvest_conc_by_testing_day.csv", header = TRUE)
shelf_life <- read.csv("data/wrangled/shelf_life_conc_by_testing_day.csv", header = TRUE)
weather_met <- read.csv("data/wrangled/weather_metrics.csv", header = TRUE)

## ----------------------Wrangling----------------------------------------------

#Defining primary growth models 

#Baranyi without lag
baranyi_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10nmax - log10(1 + (10^(log10nmax - log10n0) - 1) * exp(-mumax * day))
}

#Buchanan without lag
buchanan_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10n0 + (day <= ((log10nmax - log10n0) * log(10) / mumax)) * mumax * day / log(10) + (day > ((log10nmax - log10n0) * log(10) / mumax)) * (log10nmax - log10n0)
}

#Wrapper, for growth models

log10N <- function(day, log10n0, log10nmax, mumax, model_name) {
  if (model_name == "baranyi_no_lag") {
    return(baranyi_no_lag(day, log10n0, log10nmax, mumax))
  }
  else if(model_name == "buchanan_no_lag") {
    return(buchanan_no_lag(day, log10n0, log10nmax, mumax))
  }else {
    stop(paste0(model_name, " is not a valid model name. Must be one of baranyi_no_lag or buchanan_no_lag"))
  }
}

#Assigning locations to the growth parameters dataframe
shelf_life_location <- filter(metadata, sample_type == "f")
  
param$loc <- vector(mode = "logical", length = nrow(param))

for(i in 1:nrow(param)){
  index <- which(param$sampling_code[i] == shelf_life_location$sampling_code)
  param$loc[i] <- shelf_life_location$location[index]
}

#Filtering out parameters for samples from from Georgia (1122_B) and Florida (1221_A)
#Filtering out parameters from 0422_B. Unrealistically high Nmax; sample was only tested until day 14 of shelf life
#Growth model for 0522_A not present, since this sample was tested until day 7 of shelf life

param <- param %>%
  filter(!(sampling_code %in% c("0422_B", "0522_A", "0922_B")))

wilcox_test_param_apc <- param %>%
  filter(test == "apc") %>%
  filter(loc %in% c("AZ", "CA"))

wilcox_test_param_pc <- param %>%
  filter(test == "pc") %>%
  filter(loc %in% c("AZ", "CA"))

param_summary_median <- param %>%
  group_by(loc, test) %>% 
  summarize(median_n0 = round(quantile(n0, type = 7)[3], 2), 
            lb_n0 = round(quantile(n0, type = 7)[2], 2),
            ub_n0 = round(quantile(n0, type = 7)[4], 2),
            median_mumax = round(quantile(mumax, type = 7)[3], 2), 
            lb_mumax = round(quantile(mumax, type = 7)[2], 2),
            ub_mumax = round(quantile(mumax, type =7)[4], 2),
            median_nmax = round(quantile(nmax, type = 7)[3], 2), 
            lb_nmax = round(quantile(nmax, type = 7)[2], 2),
            ub_nmax = round(quantile(nmax, type = 7)[4], 2), 
            n = n())

## ----------------------APC----------------------------------------------------

#N0
#Checking for normality
hist(wilcox_test_param_apc$n0[wilcox_test_param_apc$loc == "AZ"], )
hist(wilcox_test_param_apc$n0[wilcox_test_param_apc$loc == "CA"], )

#Will use a Wilcoxon rank sum test, since the data is not normally distributed
n0_apc <- wilcox.test(wilcox_test_param_apc$n0[wilcox_test_param_apc$loc == "AZ"], 
                      wilcox_test_param_apc$n0[wilcox_test_param_apc$loc == "CA"])
n0_apc

#Mumuax

#Checking for normality
hist(wilcox_test_param_apc$mumax[wilcox_test_param_apc$loc == "AZ"], )
hist(wilcox_test_param_apc$mumax[wilcox_test_param_apc$loc == "CA"], )

#Will use a Wilcoxon rank sum test, since the data is not normally distributed
mumax_apc <- wilcox.test(wilcox_test_param_apc$mumax[wilcox_test_param_apc$loc == "AZ"], 
                         wilcox_test_param_apc$mumax[wilcox_test_param_apc$loc == "CA"])
mumax_apc

#Nmax

#Checking for normality
hist(wilcox_test_param_apc$nmax[wilcox_test_param_apc$loc == "AZ"], )
hist(wilcox_test_param_apc$nmax[wilcox_test_param_apc$loc == "CA"], )

#Will use a Wilcoxon rank sum test, since the data is not normally distributed
nmax_apc <- wilcox.test(wilcox_test_param_apc$nmax[wilcox_test_param_apc$loc == "AZ"], 
                        wilcox_test_param_apc$nmax[wilcox_test_param_apc$loc == "CA"])
nmax_apc

## ----------------------PC-----------------------------------------------------

#N0
#Checking for normality
hist(wilcox_test_param_pc$n0[wilcox_test_param_pc$loc == "AZ"], )
hist(wilcox_test_param_pc$n0[wilcox_test_param_pc$loc == "CA"], )

#Will use a Wilcoxon rank sum test, since the data is not normally distributed
n0_pc <- wilcox.test(wilcox_test_param_pc$n0[wilcox_test_param_pc$loc == "AZ"], 
                      wilcox_test_param_pc$n0[wilcox_test_param_pc$loc == "CA"])
n0_pc

#Mumuax

#Checking for normality
hist(wilcox_test_param_pc$mumax[wilcox_test_param_pc$loc == "AZ"], )
hist(wilcox_test_param_pc$mumax[wilcox_test_param_pc$loc == "CA"], )

#Will use a Wilcoxon rank sum test, since the data is not normally distributed
mumax_pc <- wilcox.test(wilcox_test_param_pc$mumax[wilcox_test_param_pc$loc == "AZ"], 
                         wilcox_test_param_pc$mumax[wilcox_test_param_apc$loc == "CA"])
mumax_pc

#Nmax

#Checking for normality
hist(wilcox_test_param_pc$nmax[wilcox_test_param_pc$loc == "AZ"], )
hist(wilcox_test_param_pc$nmax[wilcox_test_param_pc$loc == "CA"], )

#Will use a Wilcoxon rank sum test, since the data is not normally distributed
nmax_pc <- wilcox.test(wilcox_test_param_pc$nmax[wilcox_test_param_pc$loc == "AZ"], 
                        wilcox_test_param_pc$nmax[wilcox_test_param_pc$loc == "CA"])
nmax_pc

## ----------------------Multiple comparisons adjustment------------------------

#N0:

p_val_param_n0 <- c(n0_apc[[3]],
                 n0_pc[[3]])
p.adjust(p_val_param_n0, method = "holm", n = 2)


p_val_param_mumax <- c(mumax_apc[[3]],
                    mumax_pc[[3]])
p.adjust(p_val_param_mumax, method = "holm", n = 2)


p_val_param_nmax <- c(nmax_apc[[3]],
                    nmax_pc[[3]])
p.adjust(p_val_param_nmax, method = "holm", n = 2)

## ----------------------Time to 7 log------------------------------------------

param$time_to_7_log <- vector(mode = "logical", length = nrow(param))
param$time_to_7_log <- as.numeric(param$time_to_7_log )

day <- c(1:30)
conc <- vector(mode = "logical", length = length(day))

pred <- bind_cols(day, conc)
colnames(pred) <- c("day", "conc")

for(i in 1:nrow(param)){
  if(param$model[i] == "baranyi_no_lag"){
    my_model <- list(
      model ="Baranyi",
      logN0 = param$n0[i],
      logNmax = param$nmax[i],
      mu = param$mumax[i]/2.303,
      lambda = 0
    )
  my_pred <- predict_growth(day, 
                            my_model)
  } else {
    my_model = list(
      model ="Trilinear",
      logN0 = param$n0[i],
      logNmax = param$nmax[i],
      mu = param$mumax[i]/2.303,
      lambda = 0)
      
    my_pred <- predict_growth(day, 
                              my_model)
  }
  
  time_to_7 <- time_to_size(my_pred, 7) 
  param$time_to_7_log[i] <- time_to_7
}

#Time to 7 log 
param %>%
  group_by(loc, test) %>%
  filter(!is.na(time_to_7_log)) %>%
  summarize(median = quantile(time_to_7_log, type = 7)[3], 
            lb = quantile(time_to_7_log, type = 7)[2],
            ub = quantile(time_to_7_log, type = 7)[4])

param_apc_time <- param %>%
  filter(test == "apc") %>%
  filter(loc %in% c("AZ", "CA"))

param_pc_time <- param %>%
  filter(test == "pc") %>%
  filter(loc %in% c("AZ", "CA"))

hist(param_apc_time$time_to_7_log)

time_apc <- wilcox.test(param_apc_time$time_to_7_log[param_apc_time$loc == "AZ"], 
                        param_apc_time$time_to_7_log[param_apc_time$loc == "CA"])

time_pc <- wilcox.test(param_pc_time$time_to_7_log[param_pc_time$loc == "AZ"], 
                        param_pc_time$time_to_7_log[param_pc_time$loc == "CA"])

p_val_time_to_7_log <- c(time_apc[[3]],
                         time_pc[[3]])
p.adjust(p_val_time_to_7_log, method = "holm", n = 2)

## ----------------------Time to 7 log, standard N0-----------------------------

param$time_to_7_log_st <- vector(mode = "logical", length = nrow(param))
param$time_to_7_log_st <- as.numeric(param$time_to_7_log )

day <- c(1:30)
conc <- vector(mode = "logical", length = length(day))

pred <- bind_cols(day, conc)
colnames(pred) <- c("day", "conc")

for(i in 1:nrow(param)){
  if(param$model[i] == "baranyi_no_lag"){
    my_model <- list(
      model ="Baranyi",
      logN0 = 5,
      logNmax = param$nmax[i],
      mu = param$mumax[i]/2.303,
      lambda = 0
    )
    my_pred <- predict_growth(day, 
                              my_model)
  } else {
    my_model = list(
      model ="Trilinear",
      logN0 = 5,
      logNmax = param$nmax[i],
      mu = param$mumax[i]/2.303,
      lambda = 0)
    
    my_pred <- predict_growth(day, 
                              my_model)
  }
  
  time_to_7 <- time_to_size(my_pred, 7) 
  param$time_to_7_log_st[i] <- time_to_7
}

#Time to 7 log 
param %>%
  group_by(loc, test) %>%
  filter(!is.na(time_to_7_log_st)) %>%
  summarize(median = quantile(time_to_7_log_st, type = 7)[3], 
            lb = quantile(time_to_7_log_st, type = 7)[2],
            ub = quantile(time_to_7_log_st, type = 7)[4])

param_apc_time <- param %>%
  filter(test == "apc") %>%
  filter(loc %in% c("AZ", "CA"))

param_pc_time <- param %>%
  filter(test == "pc") %>%
  filter(loc %in% c("AZ", "CA"))

hist(param_apc_time$time_to_7_log_st)

time_st_apc <- wilcox.test(param_apc_time$time_to_7_log_st[param_apc_time$loc == "AZ"], 
                        param_apc_time$time_to_7_log_st[param_apc_time$loc == "CA"])

time_st_pc <- wilcox.test(param_pc_time$time_to_7_log_st[param_pc_time$loc == "AZ"], 
                       param_pc_time$time_to_7_log_st[param_pc_time$loc == "CA"])

p_val_time_to_7_log_st <- c(time_st_apc[[3]],
                         time_st_pc[[3]])
p.adjust(p_val_time_to_7_log_st, method = "holm", n = 2)

## ----------------------Correlation of N0 and Nmax-----------------------------


apc_n0_nmax_cor <- cor.test(param$n0[param$test == "apc"], 
         param$nmax[param$test == "apc"],
         method = "spearman")

pc_n0_nmax_cor <- cor.test(param$n0[param$test == "pc"], 
         param$nmax[param$test == "pc"],
         method = "spearman")

p_val_time_n0_nmax <- c(apc_n0_nmax_cor[[3]],
                        pc_n0_nmax_cor[[3]])

p.adjust(p_val_time_n0_nmax, method = "holm", n = 2)

## ----------------------Impact of weather on Nmax------------------------------

param_weat <- param %>%
  dplyr::select(sampling_code, nmax, test, loc)

cluster <- vector(mode = "numeric", length = 26)
var_name <- vector(mode = "character", length = 26)
aic <- vector(mode = "character", length = 26)

finish_weather <- bind_cols(cluster, var_name, aic, aic)
colnames(finish_weather) <- c("cluster", "var_name", "aic_apc", "aic_pc")

finish_weather$cluster <- c(rep(1, times = 6),
                            2:4,
                            rep(5, times = 8),
                            rep(6, times = 9))

finish_weather$var_name <- c("windspeed_mean_1d",  "windspeed_mean_2d",
                             "srad_mean_1d", "srad_mean_2d", 
                             "srad_mean_3d", "srad_mean_4d_7d",
                             "precip_sum_1d_3d", "windspeed_mean_4d_7d",
                             "windspeed_mean_3d", "temp_mean_1d",
                             "temp_min_1d", "temp_max_1d",
                             "dew_mean_1d", "dew_mean_2d",
                             "dew_mean_3d", "dew_mean_4d_7d",
                             "precip_sum_4d_7d", "temp_mean_2d",
                             "temp_mean_3d", "temp_mean_4d_7d",
                             "temp_min_2d", "temp_min_3d", 
                             "temp_min_4d_7d", "temp_max_2d",
                             "temp_max_3d", "temp_max_4d_7d")

finished_met <- weather_met %>%
  filter(sample_type != "r")

param_weather <- merge(param_weat,
                                     finished_met,
                                     by.x = "sampling_code",
                                     by.y = "sampling_code")

param_weather_apc <- filter(param_weather, test == "apc")
param_weather_pc <- filter(param_weather, test == "pc")



for(i in 1:nrow(finish_weather)){
  form <-  formula(paste("nmax ~ loc + ", finish_weather$var_name[i]))
  f_lin_mod_weath_apc <- lm(form, 
                            data = filter(param_weather_apc, test == "apc"))
  f_lin_mod_weath_pc <-lm(form, 
                          data = filter(param_weather_pc, test == "pc"))
  finish_weather$aic_apc[i] <- AIC(f_lin_mod_weath_apc)
  finish_weather$aic_pc[i] <- AIC(f_lin_mod_weath_pc)
}

for(i in 1:nrow(finish_weather)){
  form <-  formula(paste("nmax ~ loc + ", finish_weather$var_name[i]))
  f_lin_mod_weath_apc <- lm(form, 
                                        data = filter(param_weather_apc, test == "apc"))
  f_lin_mod_weath_pc <-lm(form, 
                                       data = filter(param_weather_pc, test == "pc"))
  finish_weather$aic_apc[i] <- AIC(f_lin_mod_weath_apc)
  finish_weather$aic_pc[i] <- AIC(f_lin_mod_weath_pc)
}

finish_weather$aic_apc <- as.numeric(finish_weather$aic_apc)
finish_weather$aic_pc <- as.numeric(finish_weather$aic_pc)

#APC, weather variable selection:
#APC cluster 1: srad_mean_3d, windspeed_mean_1d
cor(param_weather_apc$srad_mean_3d,
    param_weather_apc$windspeed_mean_1d, 
    method = "spearman")

#APC cluster 5: dew_mean_2d, temp_mean_1d
cor(param_weather_apc$dew_mean_2d,
    param_weather_apc$dew_mean_3d, 
    method = "spearman")

cor(param_weather_apc$dew_mean_2d,
    param_weather_apc$temp_min_1d, 
    method = "spearman")

cor(param_weather_apc$dew_mean_2d,
    param_weather_apc$temp_mean_1d, 
    method = "spearman")

#APC cluster 6: temp_min_4d_7d, temp_max_3d
cor(param_weather_apc$temp_min_4d_7d,
    param_weather_apc$temp_min_2d, 
    method = "spearman")

cor(param_weather_apc$temp_min_4d_7d,
    param_weather_apc$temp_mean_4d_7d, 
    method = "spearman")

cor(param_weather_apc$temp_min_4d_7d,
    param_weather_apc$temp_mean_3d, 
    method = "spearman")

cor(param_weather_apc$temp_min_4d_7d,
    param_weather_apc$temp_mean_2d, 
    method = "spearman")


cor(param_weather_apc$temp_min_4d_7d,
    param_weather_apc$temp_max_3d, 
    method = "spearman")

apc_weather_var <- finish_weather %>%
  filter(var_name %in% c("srad_mean_3d", "windspeed_mean_1d", "dew_mean_2d", 
                                "temp_mean_1d", "temp_min_4d_7d", "temp_max_3d",
                                "precip_sum_1d_3d", "windspeed_mean_4d_7d", 
                         "windspeed_mean_3d"))

cor_apc_dat <- param_weather_apc %>%
  select(srad_mean_3d, windspeed_mean_1d, dew_mean_2d,
         temp_mean_1d, temp_min_4d_7d,
         precip_sum_1d_3d, temp_max_3d,
         windspeed_mean_4d_7d, windspeed_mean_3d)

cor_apc <- cor(cor_apc_dat, method = "spearman")
#Replacing dew_min_2d, due to high correlation (>0.70) with temp_min_4d_7d
#dew_min_2d -> temp_max_1d (see below); could not choose other variables
# that led to a low model AIC (e.g., temp_min_1d) due to correlation (>0.70) with other weather variables

#APC cluster 5: temp_mean_1d,dew_mean_4d_7d
cor(param_weather_apc$temp_mean_1d,
    param_weather_apc$temp_min_1d, 
    method = "spearman")

cor(param_weather_apc$dew_mean_1d,
    param_weather_apc$temp_mean_1d, 
    method = "spearman")

cor(param_weather_apc$temp_mean_1d,
    param_weather_apc$temp_max_1d, 
    method = "spearman")

cor(param_weather_apc$temp_mean_1d,
    param_weather_apc$dew_mean_4d_7d, 
    method = "spearman")


apc_weather_var <- finish_weather %>%
  filter(var_name %in% c("srad_mean_3d", "windspeed_mean_1d", "temp_mean_1d", 
                         "dew_mean_4d_7d", "temp_min_4d_7d", "temp_max_3d",
                         "precip_sum_1d_3d", "windspeed_mean_4d_7d", 
                         "windspeed_mean_3d"))

cor_apc_dat <- param_weather_apc %>%
  select(srad_mean_3d, windspeed_mean_1d, temp_mean_1d, 
         dew_mean_4d_7d, temp_min_4d_7d, temp_max_3d,
         precip_sum_1d_3d, windspeed_mean_4d_7d, 
         windspeed_mean_3d)

cor_apc <- cor(cor_apc_dat, method = "spearman")

#Order to add to model: temp_min_4d_7d, precip_sum_1d_3d, srad_mean_3d, 
#temp_max_3d, windspeed_mean_1d, temp_mean_1d, windspeed_mean_4d_7d, windspeed_mean_3d, 
#dew_mean_4d_7d

nmax_mod_1_apc <- lm(nmax ~ loc + temp_min_4d_7d, data = param_weather_apc)

step(nmax_mod_1_apc, test = "F")
summary(nmax_mod_1_apc)

#PC, weather variable selection: 
#PC cluster 1: srad_mean_3d, windspeed_mean_1d
cor(param_weather_pc$srad_mean_3d,
    param_weather_pc$windspeed_mean_1d, 
    method = "spearman")

#PC cluster 5: temp_min_1d, temp_max_1d
cor(param_weather_pc$temp_min_1d,
    param_weather_pc$dew_mean_2d, 
    method = "spearman")

cor(param_weather_pc$temp_min_1d,
    param_weather_pc$temp_mean_1d, 
    method = "spearman")

cor(param_weather_pc$temp_min_1d,
    param_weather_pc$dew_mean_3d, 
    method = "spearman")

cor(param_weather_pc$temp_min_1d,
    param_weather_pc$temp_max_1d, 
    method = "spearman")

#PC cluster 6: temp_min_4d_7d, temp_max_3d
cor(param_weather_pc$temp_min_4d_7d,
    param_weather_pc$temp_min_2d, 
    method = "spearman")

cor(param_weather_pc$temp_min_4d_7d,
    param_weather_pc$temp_max_3d, 
    method = "spearman")

pc_weather_var <- finish_weather %>%
  filter(var_name %in% c("srad_mean_3d", "windspeed_mean_1d", "temp_min_1d", 
                         "temp_max_1d", "temp_min_4d_7d", "temp_max_3d",
                         "precip_sum_1d_3d", "windspeed_mean_4d_7d", 
                         "windspeed_mean_3d"))

cor_pc_dat <- param_weather_pc %>%
  select(srad_mean_3d, windspeed_mean_1d, temp_min_1d, 
         temp_max_1d, temp_min_4d_7d, temp_max_3d,
         precip_sum_1d_3d, windspeed_mean_4d_7d, 
         windspeed_mean_3d)

cor_pc <- cor(cor_pc_dat, method = "spearman")

#Need to replace temp_min_1d due to high correlation (>0.70);

#PC cluster 5: dew_mean_3d, temp_max_1d
cor(param_weather_pc$dew_mean_2d, #dew_mean_2d is correlated with variables from other clusters
    param_weather_pc$temp_max_1d, 
    method = "spearman")

cor(param_weather_pc$dew_mean_3d, 
    param_weather_pc$temp_max_1d, 
    method = "spearman")

pc_weather_var <- finish_weather %>%
  filter(var_name %in% c("srad_mean_3d", "windspeed_mean_1d", "dew_mean_3d", 
                         "temp_max_1d", "temp_min_4d_7d", "temp_max_3d",
                         "precip_sum_1d_3d", "windspeed_mean_4d_7d", 
                         "windspeed_mean_3d"))

cor_pc_dat <- param_weather_pc %>%
  select(srad_mean_3d, windspeed_mean_1d, dew_mean_3d, 
         temp_max_1d, temp_min_4d_7d, temp_max_3d,
         precip_sum_1d_3d, windspeed_mean_4d_7d, 
         windspeed_mean_3d)

cor_pc <- cor(cor_pc_dat, method = "spearman")

#Order of adding variables:temp_min_4d_7d, precip_sum_1d_3d, srad_mean_3d, 
#temp_max_3d, windspeed_mean_1d, dew_mean_3d, windspeed_mean_4d_7d, 
#temp_max_1d, windspeed_mean_3d

nmax_mod_1_pc <- lm(nmax ~ loc + temp_min_4d_7d , data = param_weather_pc)

step(nmax_mod_1_pc, test = "F")
summary(nmax_mod_1_pc)

