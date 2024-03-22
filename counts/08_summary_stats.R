## ----------------------Title--------------------------------------------------
#   Statistical analysis

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

#  Script description: Calculating summary stats (e.g., median) for bacterial concentration 
#on baby spinach

## ----------------------Packages-----------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(lme4); library(ggplot2); library(lmerTest); 
library(WRS2); library(emmeans)

## ----------------------Reading in Data----------------------------------------
shelf_life_geom <- read.csv("data/wrangled/shelf_life_conc_by_testing_day.csv", header = TRUE)
harvest_geom <- read.csv("data/wrangled/harvest_conc_by_testing_day.csv", header = TRUE)
metadata <- read.csv("data/wrangled/wrangled_metadata.csv", header = TRUE)

## ----------------------Harvest overall----------------------------------------
#Filtering out data from 0722_B, and 1222_A
harvest_geom_filt <- harvest_geom %>%
  filter(!(sampling_code %in% c("0722_B"))) %>%
  filter(!(sampling_code %in% c("1222_A") & test %in% c("APC", "GN")))

harvest_geom_med_iqr_by_test <- harvest_geom_filt %>%
  group_by(test) %>%
  summarize(median = round(quantile(log_conc_geom, type = 7)[3], 2), 
            lb = round(quantile(log_conc_geom, type = 7)[2], 2), 
            ub = round(quantile(log_conc_geom, type = 7)[4], 2), 
            n = n())

harvest_geom_med_iqr_by_loc <- harvest_geom_filt %>%
  group_by(loc, test) %>%
  summarize(median = round(quantile(log_conc_geom, type = 7)[3], 2), 
            lb = round(quantile(log_conc_geom, type = 7)[2], 2), 
            ub = round(quantile(log_conc_geom, type = 7)[4], 2), 
            n = n())


## ----------------------Finished product overall-------------------------------

#Filtering out the 1122_B shelf life counts, since these samples were from GA
#Also removing counts from 0922_B (which did not have any metadata)
#Filtering out counts from 0422_B and 0522_A, because we had incomplete data for these lots
#Keeping 0422_A. This sample was tested every 7 days, until day 21 of shelf 
shelf_life_geom_filt <- shelf_life_geom %>%
  filter(!(sampling_code %in% c("0922_B")))

net_growth <- shelf_life_geom_filt %>%
  filter(!(sampling_code %in% c("0422_B", "0522_A"))) %>%
  group_by(loc,sampling_code, test) %>%
  arrange(log_conc_geom) %>%
  summarize(net_growth = log_conc_geom[n()] - log_conc_geom[1]) %>%
  ungroup()

net_growth %>%
  group_by(loc, test) %>%
  summarize(median = quantile(net_growth)[3],
            lb = quantile(net_growth)[2],
            ub = quantile(net_growth)[4])

shelf_life_geom_filt_overall <- shelf_life_geom_filt %>%
  group_by(test) %>%
  filter(day == "7") %>%
  summarize(median = round(quantile(log_conc_geom, type = 7)[3], 2), 
            lb = round(quantile(log_conc_geom, type = 7)[2], 2), 
            ub = round(quantile(log_conc_geom, type = 7)[4], 2), 
            n = n())

shelf_life_geom_med_iqr_by_loc <- shelf_life_geom_filt %>%
  group_by(loc, test) %>%
  filter(day == "7") %>%
  summarize(median = round(quantile(log_conc_geom, type = 7)[3], 2), 
            lb = round(quantile(log_conc_geom, type = 7)[2], 2), 
            ub = round(quantile(log_conc_geom, type = 7)[4], 2), 
            n = n())

## ----------------------Other--------------------------------------------------

paired_samples <- metadata %>%
  filter(paired == "yes") %>%
  group_by(sampling_code) %>%
  slice(1) %>%
  select(sampling_code) %>%
  c()

harvest_paired <- harvest_geom %>%
  filter(sampling_code %in% paired_samples$sampling_code) %>%
  select(sampling_code, test, log_conc_geom)

colnames(harvest_paired)[3] <-  "log_conc_geom_h"

d7_paired <- shelf_life_geom %>%
  filter(sampling_code %in% paired_samples$sampling_code) %>%
  filter(day == 7) %>%
  select(sampling_code, test, log_conc_geom)

colnames(d7_paired)[3] <-  "log_conc_geom_d7"

paired_data_h_d7 <- merge(harvest_paired, d7_paired,
                          by.x = c("sampling_code", "test"),
                          by.y = c("sampling_code", "test"))

hist(paired_data_h_d7$log_conc_geom_h[paired_data_h_d7$test == "APC"])
hist(paired_data_h_d7$log_conc_geom_d7[paired_data_h_d7$test == "APC"])

apc_cor <- cor.test(paired_data_h_d7$log_conc_geom_h[paired_data_h_d7$test == "APC"],
               paired_data_h_d7$log_conc_geom_d7[paired_data_h_d7$test == "APC"],
               method = "spearman")
apc_cor

pc_cor <- cor.test(paired_data_h_d7$log_conc_geom_h[paired_data_h_d7$test == "PC"],
               paired_data_h_d7$log_conc_geom_d7[paired_data_h_d7$test == "PC"],
              method = "spearman")
pc_cor

p.adjust(c(apc_cor$p.value, pc_cor$p.value), method = "holm", n = 2)

