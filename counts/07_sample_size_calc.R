## ----------------------Title--------------------------------------------------
#   Sample size calculation

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

#  Script description: Conduct a post-hoc sample size calculation using the 
#APC and PC data from the packaged samples

## ----------------------Packages-----------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(lme4); library(ggplot2); library(lmerTest); 
library(brms); library(mgcv); library(plotly)

## ----------------------Reading in Data----------------------------------------

shelf_life_geom <- read.csv("data/wrangled/shelf_life_conc_by_testing_day.csv", header = TRUE)

location_all <- read.csv("data/raw/sample_location_of_origin.csv", header = TRUE)

## ----------------------Finished product overall-------------------------------

shelf_life_geom$days_harvest_arrival <- (shelf_life_geom$hr_harvest_arrival)/24

shelf_life_geom %>%
  group_by(loc, day) %>%
  summarize(var(log_conc_geom))

shelf_life_geom_apc_samp_num <- shelf_life_geom %>%
  group_by(test, loc, sampling_code) %>%
  slice(1)

#Filtering out the 1122_B shelf life counts, since these samples were from GA
#Also removing counts from Florida (1221_A)
shelf_life_geom_apc <- shelf_life_geom %>%
  filter(test == "APC" & !(sampling_code %in% c("1122_B", "0922_B") & loc != "FL")) 

shelf_life_geom_apc$loc <- as.factor(shelf_life_geom_apc$loc)
shelf_life_geom_apc$sampling_code <- as.factor(shelf_life_geom_apc$sampling_code)
shelf_life_geom_apc$log_conc_geom <- as.numeric(shelf_life_geom_apc$log_conc_geom)
shelf_life_geom_apc$day <- as.numeric(shelf_life_geom_apc$day)

shelf_life_geom_pc <- shelf_life_geom %>%
  filter(test == "PC" & !(sampling_code %in% c("1122_B", "0922_B") & loc != "FL")) 

shelf_life_geom_pc$loc <- as.factor(shelf_life_geom_pc$loc)
shelf_life_geom_pc$sampling_code <- as.factor(shelf_life_geom_pc$sampling_code)
shelf_life_geom_pc$log_conc_geom <- as.numeric(shelf_life_geom_pc$log_conc_geom)
shelf_life_geom_pc$day <- as.numeric(shelf_life_geom_pc$day)

#Setting up a dataframe for sample size calculations
df <- crossing(ca = 4:12, az = 4:8)
df <- df %>%
  slice(rep(1:n(), each = 100))
df$iteration <- 1:nrow(df)
df$coef <- NA
df$p_val <- NA

df_apc <- df
df_pc <- df

#Setting up the sample size calculation
ca_samples <- unique(shelf_life_geom_apc$sampling_code[shelf_life_geom_apc$loc == "CA"])
az_samples <- unique(shelf_life_geom_apc$sampling_code[shelf_life_geom_apc$loc == "AZ"])

set.seed(1)
for(i in 1:nrow(df_apc)){
  ca <- sample(ca_samples, size = df_apc$ca[i], replace = FALSE)
  az <- sample(az_samples, size = df_apc$az[i], replace = FALSE)
  samps <- c(ca, az)
  
  dat <- shelf_life_geom_apc %>%
    filter(sampling_code %in% samps)
  
  mod <- lmer(log_conc_geom ~ loc + poly(day,2) + 
                 (1|sampling_code) + 
                 days_harvest_arrival, 
               REML = TRUE, data = dat)
  if(isSingular(mod) == "TRUE"){
    
    df_apc$coef[i] <- NA
    df_apc$p_val[i] <- NA
    
  } else {
    result <- summary(mod)
    
    df_apc$coef[i] <- result$coefficients[[2,2]]
    df_apc$p_val[i] <- result$coefficients[[2,5]]
  }

  
  rm(mod)
}

set.seed(1)
for(i in 1:nrow(df_pc)){
  ca <- sample(ca_samples, size = df_pc$ca[i], replace = FALSE)
  az <- sample(az_samples, size = df_pc$az[i], replace = FALSE)
  samps <- c(ca, az)
  
  dat <- shelf_life_geom_pc %>%
    filter(sampling_code %in% samps)
  
  mod <- lmer(log_conc_geom ~ loc + poly(day,2) + 
                (1|sampling_code) + 
                days_harvest_arrival, 
              REML = TRUE, data = dat)
  if(isSingular(mod) == "TRUE"){
    
    df_pc$coef[i] <- NA
    df_pc$p_val[i] <- NA
    
  } else {
    result <- summary(mod)
    
    df_pc$coef[i] <- result$coefficients[[2,2]]
    df_pc$p_val[i] <- result$coefficients[[2,5]]
  }
  
  
  rm(mod)
}


## ----------------------Summarizing--------------------------------------------

summary_apc <- df_apc %>%
  group_by(ca, az) %>%
  filter(!is.na(coef)) %>%
  summarize(mean_param = mean(coef), sd_param = sd(coef), 
            mean_p_val = mean(p_val), sd_p_val = sd(p_val),
            no_null_reject = sum(p_val < 0.05))

summary_pc <- df_pc %>%
  group_by(ca, az) %>%
  filter(!is.na(coef)) %>%
  summarize(mean_param = mean(coef), sd_param = sd(coef), 
            mean_p_val = mean(p_val), sd_p_val = sd(p_val),
            no_null_reject = sum(p_val < 0.05))
