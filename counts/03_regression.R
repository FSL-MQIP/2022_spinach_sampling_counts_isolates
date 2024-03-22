## ----------------------Title--------------------------------------------------
#   Statistical analysis

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

#  Script description: Analysis to assess whether microbial counts in Arizona
# or California are significantly different at harvest and over shelf life

#Null hypothesis: The microbial counts at harvest and throughout shelf life are 
#                 identical for spinach from Arizona and California

#Alternative hypothesis: The microbial counts at harvest and/or throughout shelf
#                 life are significantly different for spinach from Arizona and 
#                 California 

## ----------------------Packages-----------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(lme4); library(ggplot2); library(lmerTest); 
library(WRS2); library(emmeans)

## ----------------------Reading in Data----------------------------------------
shelf_life_geom <- read.csv("data/wrangled/shelf_life_conc_by_testing_day.csv", header = TRUE)
harvest_geom <- read.csv("data/wrangled/harvest_conc_by_testing_day.csv", header = TRUE)
weather_met <- read.csv("data/wrangled/weather_metrics.csv", header = TRUE)
metadata <- read.csv("data/wrangled/wrangled_metadata.csv", header = TRUE)

## ----------------------Harvest overall----------------------------------------
#Filtering out data from 1221_A, since they were from Florida; 
#1222_A was removed, because it only had data from PC
#Filtering out 0722_B, since it was delayed in shipment to Ithaca
harvest_geom_lin <- harvest_geom %>%
  filter(!(sampling_code %in% c("1221_A", "1222_A", "0722_B")))

harvest_geom_lin$loc <- as.factor(harvest_geom_lin$loc)
harvest_geom_lin$sampling_code <- as.factor(harvest_geom_lin$sampling_code)
harvest_geom_lin$log_conc_geom <- as.numeric(harvest_geom_lin$log_conc_geom)
harvest_geom_lin$test <- as.factor(harvest_geom_lin$test)

#Linear model
h_lin_mod_1 <- lmer(log_conc_geom ~ loc*test  + hr_harvest_testing + 
                      days_since_last_irrigation + (1|sampling_code), data = filter(harvest_geom_lin))
step(h_lin_mod_1) 
summary(h_lin_mod_1)

#Indicates that the following terms can be dropped: loc:test, days_since_last_irrigation

#Linear model, without weather
h_lin_mod_2 <-  lm(log_conc_geom ~ loc  + hr_harvest_testing, data = filter(harvest_geom_lin, test == "APC"))

lmerTest::step(h_lin_mod_2, test = "F")
summary(h_lin_mod_2)

#Assessing individual weather variables:

harvest_met <- weather_met %>%
  filter(sample_type != "f")

cluster <- vector(mode = "numeric", length = 26)
var_name <- vector(mode = "character", length = 26)
aic <- vector(mode = "character", length = 26)

har_weather <- bind_cols(cluster, var_name, aic)
colnames(har_weather) <- c("cluster", "var_name", "aic")

har_weather$cluster <- c(rep(1, times = 5),
                         rep(2, times = 21))

har_weather$var_name <- c("precip_sum_1d_3d","precip_sum_4d_7d",
                          "windspeed_mean_1d", "windspeed_mean_2d",
                          "windspeed_mean_3d","temp_min_2d", 
                          "temp_min_4d_7d", "temp_min_3d",
                          "temp_min_1d", "dew_mean_1d",
                          "dew_mean_2d", "dew_mean_3d",
                          "dew_mean_4d_7d", "temp_mean_3d",
                          "temp_mean_4d_7d", "temp_max_3d",
                          "temp_max_4d_7d", "temp_mean_1d",
                          "temp_mean_2d", "temp_max_1d",
                          "temp_max_2d", "windspeed_mean_4d_7d",
                          "srad_mean_1d", "srad_mean_2d", "srad_mean_3d",
                          "srad_mean_4d_7d")

harvest_geom_lin_weather <- merge(harvest_geom_lin,
                                  harvest_met,
                                  by.x = "sampling_code",
                                  by.y = "sampling_code")

for(i in 1:nrow(har_weather)){
  form <-  formula(paste("log_conc_geom ~ loc  + hr_harvest_testing + ", har_weather$var_name[i]))
  h_lin_mod_weath <- lm(form, 
                        data = filter(harvest_geom_lin_weather, test == "APC"))
  har_weather$aic[i] <- AIC(h_lin_mod_weath)
}

har_weather$aic <- as.numeric(har_weather$aic)

#Variable selection: cluster 1, precip_sum_1d_3d, and windspeed_mean_2d
cor(harvest_geom_lin_weather$precip_sum_1d_3d, 
    harvest_geom_lin_weather$windspeed_mean_2d,
    method = "spearman")

#Variable selection: cluster 2

cor(harvest_geom_lin_weather$temp_mean_1d, 
    harvest_geom_lin_weather$temp_mean_2d,
    method = "spearman")

cor(harvest_geom_lin_weather$temp_mean_1d, 
    harvest_geom_lin_weather$temp_max_1d,
    method = "spearman")

cor(harvest_geom_lin_weather$temp_mean_1d, 
    harvest_geom_lin_weather$temp_min_2d,
    method = "spearman")

#Choose temp_mean_1d, and temp_min_2d

var_harvest_weather <- har_weather %>%
  filter(var_name %in% c("temp_mean_1d", "temp_min_2d", "precip_sum_1d_3d",
                         "windspeed_mean_2d"))

cor_harv_dat <- harvest_geom_lin_weather %>%
  select(temp_mean_1d, temp_min_2d, precip_sum_1d_3d,
         windspeed_mean_2d)

cor_harv <- cor(cor_harv_dat, method = "spearman")


#Start with adding, based on the lowest AIC: temp_mean_1d, temp_min_2d, precip_sum_1d_3d, windspeed_mean_2d

h_lin_mod_4 <-  lm(log_conc_geom ~ loc + hr_harvest_testing + temp_mean_1d,
                   data = filter(harvest_geom_lin_weather, test == "APC"))

step(h_lin_mod_4, test = "F")
summary(h_lin_mod_4) 

## ----------------------Finished product overall-------------------------------

#Filtering out the 1122_B shelf life counts, since these samples were from GA
#Also removing counts from Florida (1221_A) and 0922_B (which did not have any metadata)
shelf_life_geom_lin <- shelf_life_geom %>%
  filter(!(sampling_code %in% c("1221_A", "0922_B", "1122_B")))

shelf_life_geom_lin$loc <- as.factor(shelf_life_geom_lin$loc)
shelf_life_geom_lin$sampling_code <- as.factor(shelf_life_geom_lin$sampling_code)
shelf_life_geom_lin$log_conc_geom <- as.numeric(shelf_life_geom_lin$log_conc_geom)
shelf_life_geom_lin$day <- as.numeric(shelf_life_geom_lin$day)
shelf_life_geom_lin$test <- as.factor(shelf_life_geom_lin$test)
shelf_life_geom_lin$sampling_code_samp <- paste0(shelf_life_geom_lin$sampling_code, "-", shelf_life_geom_lin$day)
shelf_life_geom_lin$sampling_code_samp <- as.factor(shelf_life_geom_lin$sampling_code_samp)

shelf_life_geom_lin$days_harvest_packaging <- (shelf_life_geom_lin$hr_harvest_packaging)/24

shelf_life_geom_lin$days_harvest_arrival <- (shelf_life_geom_lin$hr_harvest_arrival)/24

shelf_life_geom_lin$days_arrival_packaging <- (shelf_life_geom_lin$hr_arrival_packaging)/24

#Test correlation of variables 

#Time from harvest to packaging and time from harvest to arrival at processing facility. Non-linear, used "Spearman" rank correlation
cor(shelf_life_geom_lin$days_harvest_packaging, shelf_life_geom_lin$days_harvest_arrival, method = "spearman")
cor(shelf_life_geom_lin$days_arrival_packaging, shelf_life_geom_lin$days_harvest_arrival, method = "spearman")
cor(shelf_life_geom_lin$days_arrival_packaging, shelf_life_geom_lin$days_harvest_packaging, method = "spearman")
#Will use "days harvest to arrival" and "arrival to packaging" in the model, as these two had the weakest correlation
#and capture unique aspects of the supply chain

#Linear model

#Checking whether to include days_harvest_packaging or days_arrival_packaging
sl_lin_mod_1 <- lmerTest::lmer(log_conc_geom ~ loc*test + poly(day,2) + 
                       (1|sampling_code:sampling_code_samp) + 
                       days_harvest_arrival*loc + days_arrival_packaging*loc, 
                     REML = FALSE, data = shelf_life_geom_lin)
step(sl_lin_mod_1) 
summary(sl_lin_mod_1)

#APC model without weather #(i) days from arrival until packaging and (ii) the interaction of
#location with days from harvest until arrival and days from arrival until packaging was
#excluded from the model, since they were not significant

sl_lin_mod_1_apc <- lmerTest::lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                                     days_harvest_arrival*loc + days_arrival_packaging*loc, 
                    REML = FALSE, data = filter(shelf_life_geom_lin, test == "APC"))

step(sl_lin_mod_1_apc)
summary(sl_lin_mod_1_apc) 


sl_lin_mod_2_apc <- lmerTest::lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                                     days_harvest_arrival, 
                                   REML = TRUE, data = filter(shelf_life_geom_lin, test == "APC"))

step(sl_lin_mod_2_apc)
summary(sl_lin_mod_2_apc) 

#PC model, without weather

sl_lin_mod_1_pc <- lmerTest::lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                                    days_harvest_arrival*loc,
                                  REML = FALSE, 
                                  data = filter(shelf_life_geom_lin, test == "PC"))
step(sl_lin_mod_1_pc)

sl_lin_mod_2_pc <- lmerTest::lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                                    days_harvest_arrival,
                                  REML = TRUE, 
                                  data = filter(shelf_life_geom_lin, test == "PC"))
step(sl_lin_mod_2_pc)
summary(sl_lin_mod_2_pc) 

#Assessing individual weather variables:

finished_met <- weather_met %>%
  filter(sample_type != "r")

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


shelf_life_geom_lin_weather <- merge(shelf_life_geom_lin,
                                     finished_met,
                                  by.x = "sampling_code",
                                  by.y = "sampling_code")

for(i in 1:nrow(finish_weather)){
  form <-  formula(paste("log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + days_harvest_arrival + ", finish_weather$var_name[i]))
  f_lin_mod_weath_apc <- lmerTest::lmer(form, 
                        data = filter(shelf_life_geom_lin_weather, test == "APC"))
  f_lin_mod_weath_pc <- lmerTest::lmer(form, 
                            data = filter(shelf_life_geom_lin_weather, test == "PC"))
  finish_weather$aic_apc[i] <- AIC(f_lin_mod_weath_apc)
  finish_weather$aic_pc[i] <- AIC(f_lin_mod_weath_pc)
}

finish_weather$aic_apc <- as.numeric(finish_weather$aic_apc)
finish_weather$aic_pc <- as.numeric(finish_weather$aic_pc)

#APC, weather variable selection
#APC cluster 1: windspeed_mean_1d, windspeed_mean_2d
cor(shelf_life_geom_lin_weather$windspeed_mean_2d,
    shelf_life_geom_lin_weather$windspeed_mean_1d, 
    method = "spearman")

#APC cluster 5: temp_min_1d, precip_sum_4d_7d
cor(shelf_life_geom_lin_weather$temp_min_1d,
    shelf_life_geom_lin_weather$precip_sum_4d_7d, 
    method = "spearman")


#APC cluster 6: temp_min_4d_7d, temp_max_3d
cor(shelf_life_geom_lin_weather$temp_min_4d_7d,
    shelf_life_geom_lin_weather$temp_min_2d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_min_4d_7d,
    shelf_life_geom_lin_weather$temp_mean_4d_7d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_min_4d_7d,
    shelf_life_geom_lin_weather$temp_mean_2d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_min_4d_7d,
    shelf_life_geom_lin_weather$temp_min_3d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_min_4d_7d,
    shelf_life_geom_lin_weather$temp_mean_3d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_min_4d_7d,
    shelf_life_geom_lin_weather$temp_max_3d, 
    method = "spearman")

#Cor of selected variables: APC

cor_apc_dat <- shelf_life_geom_lin_weather %>%
  select(windspeed_mean_1d, windspeed_mean_2d, precip_sum_1d_3d,
         windspeed_mean_4d_7d, windspeed_mean_3d, temp_min_1d,
         precip_sum_4d_7d, temp_min_4d_7d, temp_max_3d)

cor_apc <- cor(cor_apc_dat, method = "spearman")

apc_weather_var <- finish_weather %>%
  filter(var_name %in% c("windspeed_mean_1d", "windspeed_mean_2d", "precip_sum_1d_3d",
         "windspeed_mean_4d_7d", "windspeed_mean_3d", "temp_min_1d",
         "precip_sum_4d_7d", "temp_min_4d_7d", "temp_max_3d"))

#Replacing temp_min_1d, due to high correlation (>0.70) with temp_min_4d_7d
#APC cluster 2: temp_mean_1d, precip_sum_4d_7d
cor(shelf_life_geom_lin_weather$temp_min_1d,
    shelf_life_geom_lin_weather$precip_sum_4d_7d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_mean_1d,
    shelf_life_geom_lin_weather$precip_sum_4d_7d, 
    method = "spearman")

#Cor of selected variables: APC

cor_apc_dat <- shelf_life_geom_lin_weather %>%
  select(windspeed_mean_1d, windspeed_mean_2d, precip_sum_1d_3d,
         windspeed_mean_4d_7d, windspeed_mean_3d, temp_mean_1d,
         precip_sum_4d_7d, temp_min_4d_7d, temp_max_3d)

cor_apc <- cor(cor_apc_dat, method = "spearman")

apc_weather_var <- finish_weather %>%
  filter(var_name %in% c("windspeed_mean_1d", "windspeed_mean_2d", "precip_sum_1d_3d",
                         "windspeed_mean_4d_7d", "windspeed_mean_3d", "temp_mean_1d",
                         "precip_sum_4d_7d", "temp_min_4d_7d", "temp_max_3d"))

#APC model with weather
#Order for adding variables: temp_min_4d_7d, windspeed_mean_3d, 
#precip_sum_1d_3d, precip_sum_4d_7d, temp_mean_1d, windspeed_mean_4d_7d, 
#temp_max_3d, windspeed_mean_2d, windspeed_mean_1d


sl_lin_mod_3_apc <- lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                           days_harvest_arrival + 
                           temp_min_4d_7d + windspeed_mean_3d, 
                         REML = TRUE, data = filter(shelf_life_geom_lin_weather, test == "APC"))
step(sl_lin_mod_3_apc)

summary(sl_lin_mod_3_apc) 

#PC, weather variable selection
#PC cluster 1: windspeed_mean_1d, windspeed_mean_2d
cor(shelf_life_geom_lin_weather$windspeed_mean_2d,
    shelf_life_geom_lin_weather$windspeed_mean_1d, 
    method = "spearman")

#PC cluster 5: temp_min_1d, precip_sum_4d_7d
cor(shelf_life_geom_lin_weather$temp_min_1d,
    shelf_life_geom_lin_weather$precip_sum_4d_7d, 
    method = "spearman")

#PC cluster 6: temp_min_2d, temp_mean_3d
cor(shelf_life_geom_lin_weather$temp_min_4d_7d,
    shelf_life_geom_lin_weather$temp_min_2d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_min_2d,
    shelf_life_geom_lin_weather$temp_mean_2d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_min_2d,
    shelf_life_geom_lin_weather$temp_mean_3d, 
    method = "spearman")

#Cor of selected variables: PC

cor_pc_dat <- shelf_life_geom_lin_weather %>%
  select(windspeed_mean_1d, windspeed_mean_2d, precip_sum_1d_3d,
         windspeed_mean_4d_7d, windspeed_mean_3d, temp_min_1d,
         precip_sum_4d_7d, temp_min_2d, temp_mean_3d)

cor_pc <- cor(cor_pc_dat, method = "spearman")

pc_weather_var <- finish_weather %>%
  filter(var_name %in% c("temp_min_2d", "precip_sum_1d_3d", "temp_min_1d", 
                         "precip_sum_4d_7d", "temp_mean_3d", "windspeed_mean_3d", 
                         "windspeed_mean_1d", "windspeed_mean_2d", "windspeed_mean_4d_7d"))

#Replacing temp_min_1d, as it is highly correlated (>0.70) with temp_min_2d
#PC cluster 5: temp_mean_1d, precip_sum_4d_7d
cor(shelf_life_geom_lin_weather$temp_min_1d,
    shelf_life_geom_lin_weather$precip_sum_4d_7d, 
    method = "spearman")

cor(shelf_life_geom_lin_weather$temp_mean_1d,
    shelf_life_geom_lin_weather$precip_sum_4d_7d, 
    method = "spearman")

cor_pc_dat <- shelf_life_geom_lin_weather %>%
  select(windspeed_mean_1d, windspeed_mean_2d, precip_sum_1d_3d,
         windspeed_mean_4d_7d, windspeed_mean_3d, temp_mean_1d,
         precip_sum_4d_7d, temp_min_2d, temp_mean_3d)

cor_pc <- cor(cor_pc_dat, method = "spearman")

pc_weather_var <- finish_weather %>%
  filter(var_name %in% c("temp_min_2d", "precip_sum_1d_3d", "temp_mean_1d", 
                         "precip_sum_4d_7d", "temp_mean_3d", "windspeed_mean_3d", 
                         "windspeed_mean_1d", "windspeed_mean_2d", "windspeed_mean_4d_7d"))

#PC model with weather
#Order to add variables:
#temp_min_2d, precip_sum_1d_3d, precip_sum_4d_7d,
#temp_mean_3d, temp_min_1d, windspeed_mean_3d, 
#windspeed_mean_1d, windspeed_mean_2d, windspeed_mean_4d_7d, 


sl_lin_mod_3_pc <- lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                           days_harvest_arrival + 
                           temp_min_2d, 
                         REML = TRUE, data = filter(shelf_life_geom_lin_weather, test == "PC"))
step(sl_lin_mod_3_pc)
summary(sl_lin_mod_3_pc) 

## ----------------------Finished product, with harvest data--------------------

#Filtering out observations from paired samples (i.e., data at harvest for these 
#samples)

paired_samples <- metadata %>%
  filter(paired == "yes") %>%
  group_by(sampling_code) %>%
  slice(1) %>%
  select(sampling_code) %>%
  c()


#Adjusting harvest concentration, to account from time from harvest until testing
#All samples will be adjusted to the median time from harvest until testing, 
# for the paired samples, which was: 33.25000 h

h_samples <- unique(harvest_geom$sampling_code)
quantile(filter(harvest_geom, test == "APC")$hr_harvest_testing)

shelf_life_geom_lin_hd <- shelf_life_geom_lin %>%
  filter(sampling_code %in% paired_samples$sampling_code)

shelf_life_geom_lin_hd$harvest_apc <- NA
shelf_life_geom_lin_hd$harvest_pc <- NA

for(i in 1:nrow(shelf_life_geom_lin_hd)){
  
  apc_index <- which(shelf_life_geom_lin_hd$sampling_code[i] ==
                       harvest_geom$sampling_code &
                       "APC" == harvest_geom$test)
  
  pc_index <- which(shelf_life_geom_lin_hd$sampling_code[i] ==
                      harvest_geom$sampling_code &
                       "PC" == harvest_geom$test)
  
  if(shelf_life_geom_lin_hd$sampling_code[i] == "1222_A"){
    shelf_life_geom_lin_hd$harvest_apc[i] <- NA
    shelf_life_geom_lin_hd$harvest_apc_adj[i] <- NA
  } else {
    shelf_life_geom_lin_hd$harvest_apc[i] <- harvest_geom$log_conc_geom[apc_index]
    adj_conc <- ifelse(harvest_geom$hr_harvest_testing[apc_index] < 33.25,
                       (harvest_geom$log_conc_geom[apc_index]  - 0.01*(33.25 - 
                                                                         harvest_geom$hr_harvest_testing[apc_index])),
                       (harvest_geom$log_conc_geom[apc_index]  + 0.01*(harvest_geom$hr_harvest_testing[apc_index] -
                                                                         33.25)))
    shelf_life_geom_lin_hd$harvest_apc_adj[i] <- adj_conc
  }
  shelf_life_geom_lin_hd$harvest_pc[i] <- harvest_geom$log_conc_geom[pc_index]
  
}

shelf_life_geom_lin_hd_apc <- shelf_life_geom_lin_hd %>%
  filter(test == "APC") %>%
  filter(sampling_code != "1222_A")

shelf_life_geom_lin_hd_pc <- shelf_life_geom_lin_hd %>%
  filter(test == "PC")  %>%
  filter(sampling_code != "1222_A")

#APC model, without weather and with harvest data
 
sl_lin_mod_1_apc_hd <- lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                              days_harvest_arrival + harvest_apc, 
                         REML = TRUE, data = shelf_life_geom_lin_hd_apc)

step(sl_lin_mod_1_apc_hd)
summary(sl_lin_mod_1_apc_hd) 

sl_lin_mod_2_apc_hd <- lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                              days_harvest_arrival + harvest_apc_adj, 
                            REML = TRUE, data = shelf_life_geom_lin_hd_apc)


step(sl_lin_mod_2_apc_hd)
summary(sl_lin_mod_2_apc_hd) 

#PC model, without weather and with harvest data (could not adjust PC at harvest,
#since we only had the effect of time from harvest until testing for APC data - i.e.,
#from the harvest APC model)

sl_lin_mod_1_pc_hd <- lmer(log_conc_geom ~ loc + poly(day,2) + (1|sampling_code) + 
                           days_harvest_arrival + harvest_pc, 
                         REML = TRUE, data = shelf_life_geom_lin_hd_pc)

step(sl_lin_mod_1_pc_hd)
summary(sl_lin_mod_1_pc_hd) 
