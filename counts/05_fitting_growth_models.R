## ------------------------------Title------------------------------------------
# Estimating initial parameters in order to fit growth curves to wrangled data from the 2022 sampling (data collected up till 3 January 2022)

## ------------------------------Description------------------------------------
# Project: CIDA Spinach 

# Script description: Fitting growth parameters (i.e., Baranyi without lag and 
#Buchanan without lag) to the data from the packaged samples

## ------------------------------Packages---------------------------------------
# Loading packages
library(tidyverse); library(dplyr); library(minpack.lm)

## ------------------------------Raw Data---------------------------------------

conc <- read.csv("data/wrangled/shelf_life_conc_by_testing_day.csv", header = TRUE)
initial_estimates <- read.csv("data/wrangled/initial_parameter_estimates.csv", header = TRUE)

## ------------------Defining Primary Models------------------------------------
# Defining the primary growth functions. The Baranyi (without lag), Buchanan (without lag).
#Link to nlsMicrobio on Cran: https://cran.r-project.org/web/packages/nlsMicrobio/index.html. Mumax is in log (base e).

#Baranyi without lag
baranyi_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10nmax - log10(1 + (10^(log10nmax - log10n0) - 1) * exp(-mumax * day))
}

#Buchanan without lag
buchanan_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10n0 + (day <= ((log10nmax - log10n0) * log(10) / mumax)) * mumax * day / log(10) + (day > ((log10nmax - log10n0) * log(10) / mumax)) * (log10nmax - log10n0)
}

## ------------------------------Wrangling--------------------------------------
#Filtering data by sampling and shelf life

shelf_life_conc <- filter(conc, day != "H")
shelf_life_conc$day <- as.numeric(shelf_life_conc$day)

dat1221a_apc <- filter(shelf_life_conc, sampling_code == "1221_A" & test == "APC")

dat0122a_apc <- filter(shelf_life_conc, sampling_code == "0122_A" & test == "APC")
dat0122b_apc <- filter(shelf_life_conc, sampling_code == "0122_B" & test == "APC")

dat0222a_apc <- filter(shelf_life_conc, sampling_code == "0222_A" & test == "APC")
dat0222b_apc <- filter(shelf_life_conc, sampling_code == "0222_B" & test == "APC")

dat0322a_apc <- filter(shelf_life_conc, sampling_code == "0322_A" & test == "APC")
dat0322b_apc <- filter(shelf_life_conc, sampling_code == "0322_B" & test == "APC")

dat0422a_apc <- filter(shelf_life_conc, sampling_code == "0422_A" & test == "APC")
dat0422b_apc <- filter(shelf_life_conc, sampling_code == "0422_B" & test == "APC")

dat0522a_apc <- filter(shelf_life_conc, sampling_code == "0522_A" & test == "APC")
dat0522b_apc <- filter(shelf_life_conc, sampling_code == "0522_B" & test == "APC")
dat0522c_apc <- filter(shelf_life_conc, sampling_code == "0522_C" & test == "APC")

dat0622a_apc <- filter(shelf_life_conc, sampling_code == "0622_A" & test == "APC")
dat0622b_apc <- filter(shelf_life_conc, sampling_code == "0622_B" & test == "APC")

dat0722a_apc <- filter(shelf_life_conc, sampling_code == "0722_A" & test == "APC")

dat0822a_apc <- filter(shelf_life_conc, sampling_code == "0822_A" & test == "APC")
dat0822b_apc <- filter(shelf_life_conc, sampling_code == "0822_B" & test == "APC")

dat0922a_apc <- filter(shelf_life_conc, sampling_code == "0922_A" & test == "APC")
dat0922b_apc <- filter(shelf_life_conc, sampling_code == "0922_B" & test == "APC")

dat1022a_apc <- filter(shelf_life_conc, sampling_code == "1022_A" & test == "APC")
dat1022b_apc <- filter(shelf_life_conc, sampling_code == "1022_B" & test == "APC")

dat1122a_apc <- filter(shelf_life_conc, sampling_code == "1122_A" & test == "APC")
dat1122b_apc <- filter(shelf_life_conc, sampling_code == "1122_B" & test == "APC")

dat1222a_apc <- filter(shelf_life_conc, sampling_code == "1222_A" & test == "APC")
dat1222b_apc <- filter(shelf_life_conc, sampling_code == "1222_B" & test == "APC")

dat1221a_pc <- filter(shelf_life_conc, sampling_code == "1221_A" & test == "PC")

dat0122a_pc <- filter(shelf_life_conc, sampling_code == "0122_A" & test == "PC")
dat0122b_pc <- filter(shelf_life_conc, sampling_code == "0122_B" & test == "PC")

dat0222a_pc <- filter(shelf_life_conc, sampling_code == "0222_A" & test == "PC")
dat0222b_pc <- filter(shelf_life_conc, sampling_code == "0222_B" & test == "PC")

dat0322a_pc <- filter(shelf_life_conc, sampling_code == "0322_A" & test == "PC")
dat0322b_pc <- filter(shelf_life_conc, sampling_code == "0322_B" & test == "PC")

dat0422a_pc <- filter(shelf_life_conc, sampling_code == "0422_A" & test == "PC")
dat0422b_pc <- filter(shelf_life_conc, sampling_code == "0422_B" & test == "PC")

dat0522a_pc <- filter(shelf_life_conc, sampling_code == "0522_A" & test == "PC")
dat0522b_pc <- filter(shelf_life_conc, sampling_code == "0522_B" & test == "PC")
dat0522c_pc <- filter(shelf_life_conc, sampling_code == "0522_C" & test == "PC")

dat0622a_pc <- filter(shelf_life_conc, sampling_code == "0622_A" & test == "PC")
dat0622b_pc <- filter(shelf_life_conc, sampling_code == "0622_B" & test == "PC")

dat0722a_pc <- filter(shelf_life_conc, sampling_code == "0722_A" & test == "PC")

dat0822a_pc <- filter(shelf_life_conc, sampling_code == "0822_A" & test == "PC")
dat0822b_pc <- filter(shelf_life_conc, sampling_code == "0822_B" & test == "PC")

dat0922a_pc <- filter(shelf_life_conc, sampling_code == "0922_A" & test == "PC")
dat0922b_pc <- filter(shelf_life_conc, sampling_code == "0922_B" & test == "PC")

dat1022a_pc <- filter(shelf_life_conc, sampling_code == "1022_A" & test == "PC")
dat1022b_pc <- filter(shelf_life_conc, sampling_code == "1022_B" & test == "PC")

dat1122a_pc <- filter(shelf_life_conc, sampling_code == "1122_A" & test == "PC")
dat1122b_pc <- filter(shelf_life_conc, sampling_code == "1122_B" & test == "PC")

dat1222a_pc <- filter(shelf_life_conc, sampling_code == "1222_A" & test == "PC")
dat1222b_pc <- filter(shelf_life_conc, sampling_code == "1222_B" & test == "PC")

param_estimates_apc <- initial_estimates %>%
  filter(test == "apc")

param_estimates_pc <- initial_estimates %>%
  filter(test == "pc")

## -----------------------------Fitting growth models---------------------------
#Setting up a dataframe for recording growth parameters 
growth_parameters_apc <- data.frame(matrix(data = NA, nrow = 50, ncol = 9))
colnames(growth_parameters_apc) <- c("sampling_code", "n0", "mumax", "nmax", "model", "aic", "bic", "fit", "test")

growth_parameters_pc <- data.frame(matrix(data = NA, nrow = 50, ncol = 9))
colnames(growth_parameters_pc) <- c("sampling_code", "n0", "mumax", "nmax", "model", "aic", "bic", "fit", "test")

##----1221A: APC----

index_1221a_apc <- which(param_estimates_apc$sampling_code == "1221_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat1221a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                          baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                        data = dat1221a_apc, 
                                        start=list(
                                          log10n0 = param_estimates_apc$n0[index_1221a_apc], 
                                          log10nmax = param_estimates_apc$nmax[index_1221a_apc], 
                                          mumax = (param_estimates_apc$mumax[index_1221a_apc]*2.303)),
                                        lower = c(0, 0, 0))

growth_parameters_apc$n0[1] <- summary(dat1221a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[1] <- summary(dat1221a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[1] <- summary(dat1221a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[1] <- AIC(dat1221a_baranyi_no_lag_apc)
growth_parameters_apc$bic[1] <- BIC(dat1221a_baranyi_no_lag_apc)
growth_parameters_apc$fit[1] <- "yes"
growth_parameters_apc$sampling_code[1] <- "1221_A"
growth_parameters_apc$test[1] <- "apc"
growth_parameters_apc$model[1] <-"baranyi_no_lag"



#buchanan no lag
dat1221a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                           buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                         data = dat1221a_apc, 
                                         start=list(
                                           log10n0 = param_estimates_apc$n0[index_1221a_apc], 
                                           log10nmax = param_estimates_apc$nmax[index_1221a_apc], 
                                           mumax = (param_estimates_apc$mumax[index_1221a_apc]*2.303)),
                                         lower = c(0, 0, 0))

growth_parameters_apc$n0[2] <- summary(dat1221a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[2] <- summary(dat1221a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[2] <- summary(dat1221a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[2] <- AIC(dat1221a_buchanan_no_lag_apc)
growth_parameters_apc$bic[2] <- BIC(dat1221a_buchanan_no_lag_apc)
growth_parameters_apc$fit[2] <- "yes"
growth_parameters_apc$sampling_code[2] <- "1221_A"
growth_parameters_apc$test[2] <- "apc"
growth_parameters_apc$model[2] <-"buchanan_no_lag"



##----0122A: APC----

index_0122a_apc <- which(param_estimates_apc$sampling_code == "0122_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0122a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0122a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0122a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0122a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0122a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[3] <- summary(dat0122a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[3] <- summary(dat0122a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[3] <- summary(dat0122a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[3] <- AIC(dat0122a_baranyi_no_lag_apc)
growth_parameters_apc$bic[3] <- BIC(dat0122a_baranyi_no_lag_apc)
growth_parameters_apc$fit[3] <- "yes"
growth_parameters_apc$sampling_code[3] <- "0122_A"
growth_parameters_apc$test[3] <- "apc"
growth_parameters_apc$model[3] <-"baranyi_no_lag"



#buchanan no lag
dat0122a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0122a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0122a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0122a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0122a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[4] <- summary(dat0122a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[4] <- summary(dat0122a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[4] <- summary(dat0122a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[4] <- AIC(dat0122a_buchanan_no_lag_apc)
growth_parameters_apc$bic[4] <- BIC(dat0122a_buchanan_no_lag_apc)
growth_parameters_apc$fit[4] <- "yes"
growth_parameters_apc$sampling_code[4] <- "0122_A"
growth_parameters_apc$test[4] <- "apc"
growth_parameters_apc$model[4] <-"buchanan_no_lag"

##----0122B: APC----

index_0122b_apc <- which(param_estimates_apc$sampling_code == "0122_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0122b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0122b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0122b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0122b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0122b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[5] <- summary(dat0122b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[5] <- summary(dat0122b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[5] <- summary(dat0122b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[5] <- AIC(dat0122b_baranyi_no_lag_apc)
growth_parameters_apc$bic[5] <- BIC(dat0122b_baranyi_no_lag_apc)
growth_parameters_apc$fit[5] <- "yes"
growth_parameters_apc$sampling_code[5] <- "0122_B"
growth_parameters_apc$test[5] <- "apc"
growth_parameters_apc$model[5] <-"baranyi_no_lag"



#buchanan no lag
dat0122b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0122b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0122b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0122b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0122b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[6] <- summary(dat0122b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[6] <- summary(dat0122b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[6] <- summary(dat0122b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[6] <- AIC(dat0122b_buchanan_no_lag_apc)
growth_parameters_apc$bic[6] <- BIC(dat0122b_buchanan_no_lag_apc)
growth_parameters_apc$fit[6] <- "yes"
growth_parameters_apc$sampling_code[6] <- "0122_B"
growth_parameters_apc$test[6] <- "apc"
growth_parameters_apc$model[6] <-"buchanan_no_lag"

##----0222A: APC----

index_0222a_apc <- which(param_estimates_apc$sampling_code == "0222_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0222a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0222a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0222a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0222a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0222a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[7] <- summary(dat0222a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[7] <- summary(dat0222a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[7] <- summary(dat0222a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[7] <- AIC(dat0222a_baranyi_no_lag_apc)
growth_parameters_apc$bic[7] <- BIC(dat0222a_baranyi_no_lag_apc)
growth_parameters_apc$fit[7] <- "yes"
growth_parameters_apc$sampling_code[7] <- "0222_A"
growth_parameters_apc$test[7] <- "apc"
growth_parameters_apc$model[7] <-"baranyi_no_lag"



#buchanan no lag
dat0222a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0222a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0222a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0222a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0222a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[8] <- summary(dat0222a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[8] <- summary(dat0222a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[8] <- summary(dat0222a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[8] <- AIC(dat0222a_buchanan_no_lag_apc)
growth_parameters_apc$bic[8] <- BIC(dat0222a_buchanan_no_lag_apc)
growth_parameters_apc$fit[8] <- "yes"
growth_parameters_apc$sampling_code[8] <- "0222_A"
growth_parameters_apc$test[8] <- "apc"
growth_parameters_apc$model[8] <-"buchanan_no_lag"


##----0222B: APC----

index_0222b_apc <- which(param_estimates_apc$sampling_code == "0222_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0222b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0222b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0222b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0222b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0222b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[9] <- summary(dat0222b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[9] <- summary(dat0222b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[9] <- summary(dat0222b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[9] <- AIC(dat0222b_baranyi_no_lag_apc)
growth_parameters_apc$bic[9] <- BIC(dat0222b_baranyi_no_lag_apc)
growth_parameters_apc$fit[9] <- "yes"
growth_parameters_apc$sampling_code[9] <- "0222_B"
growth_parameters_apc$test[9] <- "apc"
growth_parameters_apc$model[9] <-"baranyi_no_lag"



#buchanan no lag
dat0222b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0222b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0222b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0222b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0222b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[10] <- summary(dat0222b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[10] <- summary(dat0222b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[10] <- summary(dat0222b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[10] <- AIC(dat0222b_buchanan_no_lag_apc)
growth_parameters_apc$bic[10] <- BIC(dat0222b_buchanan_no_lag_apc)
growth_parameters_apc$fit[10] <- "yes"
growth_parameters_apc$sampling_code[10] <- "0222_B"
growth_parameters_apc$test[10] <- "apc"
growth_parameters_apc$model[10] <-"buchanan_no_lag"

##----0322A: APC----

index_0322a_apc <- which(param_estimates_apc$sampling_code == "0322_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0322a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0322a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0322a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0322a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0322a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[11] <- summary(dat0322a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[11] <- summary(dat0322a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[11] <- summary(dat0322a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[11] <- AIC(dat0322a_baranyi_no_lag_apc)
growth_parameters_apc$bic[11] <- BIC(dat0322a_baranyi_no_lag_apc)
growth_parameters_apc$fit[11] <- "yes"
growth_parameters_apc$sampling_code[11] <- "0322_A"
growth_parameters_apc$test[11] <- "apc"
growth_parameters_apc$model[11] <-"baranyi_no_lag"



#buchanan no lag
dat0322a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0322a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0322a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0322a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0322a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[12] <- summary(dat0322a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[12] <- summary(dat0322a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[12] <- summary(dat0322a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[12] <- AIC(dat0322a_buchanan_no_lag_apc)
growth_parameters_apc$bic[12] <- BIC(dat0322a_buchanan_no_lag_apc)
growth_parameters_apc$fit[12] <- "yes"
growth_parameters_apc$sampling_code[12] <- "0322_A"
growth_parameters_apc$test[12] <- "apc"
growth_parameters_apc$model[12] <-"buchanan_no_lag"

##----0322B: APC----

index_0322b_apc <- which(param_estimates_apc$sampling_code == "0322_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0322b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0322b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0322b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0322b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0322b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[13] <- summary(dat0322b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[13] <- summary(dat0322b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[13] <- summary(dat0322b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[13] <- AIC(dat0322b_baranyi_no_lag_apc)
growth_parameters_apc$bic[13] <- BIC(dat0322b_baranyi_no_lag_apc)
growth_parameters_apc$fit[13] <- "yes"
growth_parameters_apc$sampling_code[13] <- "0322_B"
growth_parameters_apc$test[13] <- "apc"
growth_parameters_apc$model[13] <-"baranyi_no_lag"



#buchanan no lag
dat0322b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0322b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0322b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0322b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0322b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[14] <- summary(dat0322b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[14] <- summary(dat0322b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[14] <- summary(dat0322b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[14] <- AIC(dat0322b_buchanan_no_lag_apc)
growth_parameters_apc$bic[14] <- BIC(dat0322b_buchanan_no_lag_apc)
growth_parameters_apc$fit[14] <- "yes"
growth_parameters_apc$sampling_code[14] <- "0322_B"
growth_parameters_apc$test[14] <- "apc"
growth_parameters_apc$model[14] <-"buchanan_no_lag"

##----0422A: APC----

index_0422a_apc <- which(param_estimates_apc$sampling_code == "0422_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0422a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0422a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0422a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0422a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0422a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[15] <- summary(dat0422a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[15] <- summary(dat0422a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[15] <- summary(dat0422a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[15] <- AIC(dat0422a_baranyi_no_lag_apc)
growth_parameters_apc$bic[15] <- BIC(dat0422a_baranyi_no_lag_apc)
growth_parameters_apc$fit[15] <- "yes"
growth_parameters_apc$sampling_code[15] <- "0422_A"
growth_parameters_apc$test[15] <- "apc"
growth_parameters_apc$model[15] <-"baranyi_no_lag"



#buchanan no lag
dat0422a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0422a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0422a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0422a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0422a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[16] <- summary(dat0422a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[16] <- summary(dat0422a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[16] <- summary(dat0422a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[16] <- AIC(dat0422a_buchanan_no_lag_apc)
growth_parameters_apc$bic[16] <- BIC(dat0422a_buchanan_no_lag_apc)
growth_parameters_apc$fit[16] <- "yes"
growth_parameters_apc$sampling_code[16] <- "0422_A"
growth_parameters_apc$test[16] <- "apc"
growth_parameters_apc$model[16] <-"buchanan_no_lag"

##----0422B: APC----

index_0422b_apc <- which(param_estimates_apc$sampling_code == "0422_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0422b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0422b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0422b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0422b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0422b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[17] <- summary(dat0422b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[17] <- summary(dat0422b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[17] <- summary(dat0422b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[17] <- AIC(dat0422b_baranyi_no_lag_apc)
growth_parameters_apc$bic[17] <- BIC(dat0422b_baranyi_no_lag_apc)
growth_parameters_apc$fit[17] <- "yes"
growth_parameters_apc$sampling_code[17] <- "0422_B"
growth_parameters_apc$test[17] <- "apc"
growth_parameters_apc$model[17] <-"baranyi_no_lag"



#buchanan no lag
dat0422b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0422b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0422b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0422b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0422b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[18] <- summary(dat0422b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[18] <- summary(dat0422b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[18] <- summary(dat0422b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[18] <- AIC(dat0422b_buchanan_no_lag_apc)
growth_parameters_apc$bic[18] <- BIC(dat0422b_buchanan_no_lag_apc)
growth_parameters_apc$fit[18] <- "no"
growth_parameters_apc$sampling_code[18] <- "0422_B"
growth_parameters_apc$test[18] <- "apc"
growth_parameters_apc$model[18] <-"buchanan_no_lag"

##----0522A: APC----

index_0522a_apc <- which(param_estimates_apc$sampling_code == "0522_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0522a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0522a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0522a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0522a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0522a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[19] <- summary(dat0522a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[19] <- summary(dat0522a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[19] <- summary(dat0522a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[19] <- AIC(dat0522a_baranyi_no_lag_apc)
growth_parameters_apc$bic[19] <- BIC(dat0522a_baranyi_no_lag_apc)
growth_parameters_apc$fit[19] <- "no"
growth_parameters_apc$sampling_code[19] <- "0522_A"
growth_parameters_apc$test[19] <- "apc"
growth_parameters_apc$model[19] <-"baranyi_no_lag"



#buchanan no lag
dat0522a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0522a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0522a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0522a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0522a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[20] <- summary(dat0522a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[20] <- summary(dat0522a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[20] <- summary(dat0522a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[20] <- AIC(dat0522a_buchanan_no_lag_apc)
growth_parameters_apc$bic[20] <- BIC(dat0522a_buchanan_no_lag_apc)
growth_parameters_apc$fit[20] <- "no"
growth_parameters_apc$sampling_code[20] <- "0522_A"
growth_parameters_apc$test[20] <- "apc"
growth_parameters_apc$model[20] <-"buchanan_no_lag"

##----0522B: APC----

index_0522b_apc <- which(param_estimates_apc$sampling_code == "0522_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0522b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0522b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0522b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0522b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0522b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[21] <- summary(dat0522b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[21] <- summary(dat0522b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[21] <- summary(dat0522b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[21] <- AIC(dat0522b_baranyi_no_lag_apc)
growth_parameters_apc$bic[21] <- BIC(dat0522b_baranyi_no_lag_apc)
growth_parameters_apc$fit[21] <- "yes"
growth_parameters_apc$sampling_code[21] <- "0522_B"
growth_parameters_apc$test[21] <- "apc"
growth_parameters_apc$model[21] <-"baranyi_no_lag"



#buchanan no lag
dat0522b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0522b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0522b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0522b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0522b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[22] <- summary(dat0522b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[22] <- summary(dat0522b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[22] <- summary(dat0522b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[22] <- AIC(dat0522b_buchanan_no_lag_apc)
growth_parameters_apc$bic[22] <- BIC(dat0522b_buchanan_no_lag_apc)
growth_parameters_apc$fit[22] <- "no"
growth_parameters_apc$sampling_code[22] <- "0522_B"
growth_parameters_apc$test[22] <- "apc"
growth_parameters_apc$model[22] <-"buchanan_no_lag"

##----0522C: APC----

index_0522c_apc <- which(param_estimates_apc$sampling_code == "0522_C" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0522c_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0522c_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0522c_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0522c_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0522c_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[23] <- summary(dat0522c_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[23] <- summary(dat0522c_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[23] <- summary(dat0522c_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[23] <- AIC(dat0522c_baranyi_no_lag_apc)
growth_parameters_apc$bic[23] <- BIC(dat0522c_baranyi_no_lag_apc)
growth_parameters_apc$fit[23] <- "yes"
growth_parameters_apc$sampling_code[23] <- "0522_C"
growth_parameters_apc$test[23] <- "apc"
growth_parameters_apc$model[23] <-"baranyi_no_lag"



#buchanan no lag
dat0522c_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0522c_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0522c_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0522c_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0522c_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[24] <- summary(dat0522c_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[24] <- summary(dat0522c_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[24] <- summary(dat0522c_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[24] <- AIC(dat0522c_buchanan_no_lag_apc)
growth_parameters_apc$bic[24] <- BIC(dat0522c_buchanan_no_lag_apc)
growth_parameters_apc$fit[24] <- "yes"
growth_parameters_apc$sampling_code[24] <- "0522_C"
growth_parameters_apc$test[24] <- "apc"
growth_parameters_apc$model[24] <-"buchanan_no_lag"

##----0622A: APC----

index_0622a_apc <- which(param_estimates_apc$sampling_code == "0622_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0622a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0622a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0622a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0622a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0622a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[25] <- summary(dat0622a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[25] <- summary(dat0622a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[25] <- summary(dat0622a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[25] <- AIC(dat0622a_baranyi_no_lag_apc)
growth_parameters_apc$bic[25] <- BIC(dat0622a_baranyi_no_lag_apc)
growth_parameters_apc$fit[25] <- "yes"
growth_parameters_apc$sampling_code[25] <- "0622_A"
growth_parameters_apc$test[25] <- "apc"
growth_parameters_apc$model[25] <-"baranyi_no_lag"



#buchanan no lag
dat0622a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0622a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0622a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0622a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0622a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[26] <- summary(dat0622a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[26] <- summary(dat0622a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[26] <- summary(dat0622a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[26] <- AIC(dat0622a_buchanan_no_lag_apc)
growth_parameters_apc$bic[26] <- BIC(dat0622a_buchanan_no_lag_apc)
growth_parameters_apc$fit[26] <- "yes"
growth_parameters_apc$sampling_code[26] <- "0622_A"
growth_parameters_apc$test[26] <- "apc"
growth_parameters_apc$model[26] <-"buchanan_no_lag"

##----0622B: APC----

index_0622b_apc <- which(param_estimates_apc$sampling_code == "0622_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0622b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0622b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0622b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0622b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0622b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[27] <- summary(dat0622b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[27] <- summary(dat0622b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[27] <- summary(dat0622b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[27] <- AIC(dat0622b_baranyi_no_lag_apc)
growth_parameters_apc$bic[27] <- BIC(dat0622b_baranyi_no_lag_apc)
growth_parameters_apc$fit[27] <- "yes"
growth_parameters_apc$sampling_code[27] <- "0622_B"
growth_parameters_apc$test[27] <- "apc"
growth_parameters_apc$model[27] <-"baranyi_no_lag"



#buchanan no lag
dat0622b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0622b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0622b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0622b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0622b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[28] <- summary(dat0622b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[28] <- summary(dat0622b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[28] <- summary(dat0622b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[28] <- AIC(dat0622b_buchanan_no_lag_apc)
growth_parameters_apc$bic[28] <- BIC(dat0622b_buchanan_no_lag_apc)
growth_parameters_apc$fit[28] <- "yes"
growth_parameters_apc$sampling_code[28] <- "0622_B"
growth_parameters_apc$test[28] <- "apc"
growth_parameters_apc$model[28] <-"buchanan_no_lag"

##----0722A: APC----

index_0722a_apc <- which(param_estimates_apc$sampling_code == "0722_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0722a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0722a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0722a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0722a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0722a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[29] <- summary(dat0722a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[29] <- summary(dat0722a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[29] <- summary(dat0722a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[29] <- AIC(dat0722a_baranyi_no_lag_apc)
growth_parameters_apc$bic[29] <- BIC(dat0722a_baranyi_no_lag_apc)
growth_parameters_apc$fit[29] <- "yes"
growth_parameters_apc$sampling_code[29] <- "0722_A"
growth_parameters_apc$test[29] <- "apc"
growth_parameters_apc$model[29] <-"baranyi_no_lag"



#buchanan no lag
dat0722a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0722a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0722a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0722a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0722a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[30] <- summary(dat0722a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[30] <- summary(dat0722a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[30] <- summary(dat0722a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[30] <- AIC(dat0722a_buchanan_no_lag_apc)
growth_parameters_apc$bic[30] <- BIC(dat0722a_buchanan_no_lag_apc)
growth_parameters_apc$fit[30] <- "yes"
growth_parameters_apc$sampling_code[30] <- "0722_A"
growth_parameters_apc$test[30] <- "apc"
growth_parameters_apc$model[30] <-"buchanan_no_lag"

##----0822A: APC----

index_0822a_apc <- which(param_estimates_apc$sampling_code == "0822_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0822a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0822a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0822a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0822a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0822a_apc]*2.303)),
                                 lower = c(0, 0, 0),
                                 control = nls.lm.control(maxiter = 100))

growth_parameters_apc$n0[31] <- summary(dat0822a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[31] <- summary(dat0822a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[31] <- summary(dat0822a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[31] <- AIC(dat0822a_baranyi_no_lag_apc)
growth_parameters_apc$bic[31] <- BIC(dat0822a_baranyi_no_lag_apc)
growth_parameters_apc$fit[31] <- "yes"
growth_parameters_apc$sampling_code[31] <- "0822_A"
growth_parameters_apc$test[31] <- "apc"
growth_parameters_apc$model[31] <-"baranyi_no_lag"



#buchanan no lag
dat0822a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0822a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0822a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0822a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0822a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[32] <- summary(dat0822a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[32] <- summary(dat0822a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[32] <- summary(dat0822a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[32] <- AIC(dat0822a_buchanan_no_lag_apc)
growth_parameters_apc$bic[32] <- BIC(dat0822a_buchanan_no_lag_apc)
growth_parameters_apc$fit[32] <- "yes"
growth_parameters_apc$sampling_code[32] <- "0822_A"
growth_parameters_apc$test[32] <- "apc"
growth_parameters_apc$model[32] <-"buchanan_no_lag"

##----0822B: APC----

index_0822b_apc <- which(param_estimates_apc$sampling_code == "0822_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0822b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0822b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0822b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0822b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0822b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[33] <- summary(dat0822b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[33] <- summary(dat0822b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[33] <- summary(dat0822b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[33] <- AIC(dat0822b_baranyi_no_lag_apc)
growth_parameters_apc$bic[33] <- BIC(dat0822b_baranyi_no_lag_apc)
growth_parameters_apc$fit[33] <- "yes"
growth_parameters_apc$sampling_code[33] <- "0822_B"
growth_parameters_apc$test[33] <- "apc"
growth_parameters_apc$model[33] <-"baranyi_no_lag"



#buchanan no lag
dat0822b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0822b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0822b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0822b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0822b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[34] <- summary(dat0822b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[34] <- summary(dat0822b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[34] <- summary(dat0822b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[34] <- AIC(dat0822b_buchanan_no_lag_apc)
growth_parameters_apc$bic[34] <- BIC(dat0822b_buchanan_no_lag_apc)
growth_parameters_apc$fit[34] <- "yes"
growth_parameters_apc$sampling_code[34] <- "0822_B"
growth_parameters_apc$test[34] <- "apc"
growth_parameters_apc$model[34] <-"buchanan_no_lag"

##----0922A: APC----

index_0922a_apc <- which(param_estimates_apc$sampling_code == "0922_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0922a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0922a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0922a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0922a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0922a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[35] <- summary(dat0922a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[35] <- summary(dat0922a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[35] <- summary(dat0922a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[35] <- AIC(dat0922a_baranyi_no_lag_apc)
growth_parameters_apc$bic[35] <- BIC(dat0922a_baranyi_no_lag_apc)
growth_parameters_apc$fit[35] <- "yes"
growth_parameters_apc$sampling_code[35] <- "0922_A"
growth_parameters_apc$test[35] <- "apc"
growth_parameters_apc$model[35] <-"baranyi_no_lag"



#buchanan no lag
dat0922a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0922a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0922a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0922a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0922a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[36] <- summary(dat0922a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[36] <- summary(dat0922a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[36] <- summary(dat0922a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[36] <- AIC(dat0922a_buchanan_no_lag_apc)
growth_parameters_apc$bic[36] <- BIC(dat0922a_buchanan_no_lag_apc)
growth_parameters_apc$fit[36] <- "yes"
growth_parameters_apc$sampling_code[36] <- "0922_A"
growth_parameters_apc$test[36] <- "apc"
growth_parameters_apc$model[36] <-"buchanan_no_lag"

##----0922B: APC----

index_0922b_apc <- which(param_estimates_apc$sampling_code == "0922_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat0922b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat0922b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_0922b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_0922b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_0922b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[37] <- summary(dat0922b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[37] <- summary(dat0922b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[37] <- summary(dat0922b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[37] <- AIC(dat0922b_baranyi_no_lag_apc)
growth_parameters_apc$bic[37] <- BIC(dat0922b_baranyi_no_lag_apc)
growth_parameters_apc$fit[37] <- "yes"
growth_parameters_apc$sampling_code[37] <- "0922_B"
growth_parameters_apc$test[37] <- "apc"
growth_parameters_apc$model[37] <-"baranyi_no_lag"

#buchanan no lag
dat0922b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat0922b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_0922b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_0922b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_0922b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[38] <- summary(dat0922b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[38] <- summary(dat0922b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[38] <- summary(dat0922b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[38] <- AIC(dat0922b_buchanan_no_lag_apc)
growth_parameters_apc$bic[38] <- BIC(dat0922b_buchanan_no_lag_apc)
growth_parameters_apc$fit[38] <- "yes"
growth_parameters_apc$sampling_code[38] <- "0922_B"
growth_parameters_apc$test[38] <- "apc"
growth_parameters_apc$model[38] <-"buchanan_no_lag"

##----1022A: APC----

index_1022a_apc <- which(param_estimates_apc$sampling_code == "1022_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat1022a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat1022a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_1022a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_1022a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_1022a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[39] <- summary(dat1022a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[39] <- summary(dat1022a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[39] <- summary(dat1022a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[39] <- AIC(dat1022a_baranyi_no_lag_apc)
growth_parameters_apc$bic[39] <- BIC(dat1022a_baranyi_no_lag_apc)
growth_parameters_apc$fit[39] <- "yes"
growth_parameters_apc$sampling_code[39] <- "1022_A"
growth_parameters_apc$test[39] <- "apc"
growth_parameters_apc$model[39] <-"baranyi_no_lag"

#buchanan no lag
dat1022a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat1022a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_1022a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_1022a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_1022a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[40] <- summary(dat1022a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[40] <- summary(dat1022a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[40] <- summary(dat1022a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[40] <- AIC(dat1022a_buchanan_no_lag_apc)
growth_parameters_apc$bic[40] <- BIC(dat1022a_buchanan_no_lag_apc)
growth_parameters_apc$fit[40] <- "yes"
growth_parameters_apc$sampling_code[40] <- "1022_A"
growth_parameters_apc$test[40] <- "apc"
growth_parameters_apc$model[40] <-"buchanan_no_lag"

##----1022B: APC----

index_1022b_apc <- which(param_estimates_apc$sampling_code == "1022_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat1022b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat1022b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_1022b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_1022b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_1022b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[41] <- summary(dat1022b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[41] <- summary(dat1022b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[41] <- summary(dat1022b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[41] <- AIC(dat1022b_baranyi_no_lag_apc)
growth_parameters_apc$bic[41] <- BIC(dat1022b_baranyi_no_lag_apc)
growth_parameters_apc$fit[41] <- "yes"
growth_parameters_apc$sampling_code[41] <- "1022_B"
growth_parameters_apc$test[41] <- "apc"
growth_parameters_apc$model[41] <-"baranyi_no_lag"

#buchanan no lag
dat1022b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat1022b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_1022b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_1022b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_1022b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[42] <- summary(dat1022b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[42] <- summary(dat1022b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[42] <- summary(dat1022b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[42] <- AIC(dat1022b_buchanan_no_lag_apc)
growth_parameters_apc$bic[42] <- BIC(dat1022b_buchanan_no_lag_apc)
growth_parameters_apc$fit[42] <- "yes"
growth_parameters_apc$sampling_code[42] <- "1022_B"
growth_parameters_apc$test[42] <- "apc"
growth_parameters_apc$model[42] <-"buchanan_no_lag"

##----1122A: APC----

index_1122a_apc <- which(param_estimates_apc$sampling_code == "1122_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat1122a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat1122a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_1122a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_1122a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_1122a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[43] <- summary(dat1122a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[43] <- summary(dat1122a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[43] <- summary(dat1122a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[43] <- AIC(dat1122a_baranyi_no_lag_apc)
growth_parameters_apc$bic[43] <- BIC(dat1122a_baranyi_no_lag_apc)
growth_parameters_apc$fit[43] <- "yes"
growth_parameters_apc$sampling_code[43] <- "1122_A"
growth_parameters_apc$test[43] <- "apc"
growth_parameters_apc$model[43] <-"baranyi_no_lag"

#buchanan no lag
dat1122a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat1122a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_1122a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_1122a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_1122a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[44] <- summary(dat1122a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[44] <- summary(dat1122a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[44] <- summary(dat1122a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[44] <- AIC(dat1122a_buchanan_no_lag_apc)
growth_parameters_apc$bic[44] <- BIC(dat1122a_buchanan_no_lag_apc)
growth_parameters_apc$fit[44] <- "yes"
growth_parameters_apc$sampling_code[44] <- "1122_A"
growth_parameters_apc$test[44] <- "apc"
growth_parameters_apc$model[44] <-"buchanan_no_lag"

##----1122B: APC----

index_1122b_apc <- which(param_estimates_apc$sampling_code == "1122_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat1122b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat1122b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_1122b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_1122b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_1122b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[45] <- summary(dat1122b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[45] <- summary(dat1122b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[45] <- summary(dat1122b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[45] <- AIC(dat1122b_baranyi_no_lag_apc)
growth_parameters_apc$bic[45] <- BIC(dat1122b_baranyi_no_lag_apc)
growth_parameters_apc$fit[45] <- "yes"
growth_parameters_apc$sampling_code[45] <- "1122_B"
growth_parameters_apc$test[45] <- "apc"
growth_parameters_apc$model[45] <-"baranyi_no_lag"

#buchanan no lag
dat1122b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat1122b_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_1122b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_1122b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_1122b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[46] <- summary(dat1122b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[46] <- summary(dat1122b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[46] <- summary(dat1122b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[46] <- AIC(dat1122b_buchanan_no_lag_apc)
growth_parameters_apc$bic[46] <- BIC(dat1122b_buchanan_no_lag_apc)
growth_parameters_apc$fit[46] <- "yes"
growth_parameters_apc$sampling_code[46] <- "1122_B"
growth_parameters_apc$test[46] <- "apc"
growth_parameters_apc$model[46] <-"buchanan_no_lag"

##----1222A: APC----

index_1222a_apc <- which(param_estimates_apc$sampling_code == "1222_A" & param_estimates_apc$test == "apc")

#baranyi no lag
dat1222a_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat1222a_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_1222a_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_1222a_apc], 
                                   mumax = (param_estimates_apc$mumax[index_1222a_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[47] <- summary(dat1222a_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[47] <- summary(dat1222a_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[47] <- summary(dat1222a_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[47] <- AIC(dat1222a_baranyi_no_lag_apc)
growth_parameters_apc$bic[47] <- BIC(dat1222a_baranyi_no_lag_apc)
growth_parameters_apc$fit[47] <- "yes"
growth_parameters_apc$sampling_code[47] <- "1222_A"
growth_parameters_apc$test[47] <- "apc"
growth_parameters_apc$model[47] <-"baranyi_no_lag"

#buchanan no lag
dat1222a_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat1222a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_1222a_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_1222a_apc], 
                                    mumax = (param_estimates_apc$mumax[index_1222a_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[48] <- summary(dat1222a_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[48] <- summary(dat1222a_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[48] <- summary(dat1222a_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[48] <- AIC(dat1222a_buchanan_no_lag_apc)
growth_parameters_apc$bic[48] <- BIC(dat1222a_buchanan_no_lag_apc)
growth_parameters_apc$fit[48] <- "yes"
growth_parameters_apc$sampling_code[48] <- "1222_A"
growth_parameters_apc$test[48] <- "apc"
growth_parameters_apc$model[48] <-"buchanan_no_lag"

##----1222B: APC----

index_1222b_apc <- which(param_estimates_apc$sampling_code == "1222_B" & param_estimates_apc$test == "apc")

#baranyi no lag
dat1222b_baranyi_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                   baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                 data = dat1222b_apc, 
                                 start=list(
                                   log10n0 = param_estimates_apc$n0[index_1222b_apc], 
                                   log10nmax = param_estimates_apc$nmax[index_1222b_apc], 
                                   mumax = (param_estimates_apc$mumax[index_1222b_apc]*2.303)),
                                 lower = c(0, 0, 0))

growth_parameters_apc$n0[49] <- summary(dat1222b_baranyi_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[49] <- summary(dat1222b_baranyi_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[49] <- summary(dat1222b_baranyi_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[49] <- AIC(dat1222b_baranyi_no_lag_apc)
growth_parameters_apc$bic[49] <- BIC(dat1222b_baranyi_no_lag_apc)
growth_parameters_apc$fit[49] <- "yes"
growth_parameters_apc$sampling_code[49] <- "1222_B"
growth_parameters_apc$test[49] <- "apc"
growth_parameters_apc$model[49] <-"baranyi_no_lag"

#buchanan no lag
dat1222b_buchanan_no_lag_apc <- nlsLM(log_conc_geom ~ 
                                    buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                  data = dat1222a_apc, 
                                  start=list(
                                    log10n0 = param_estimates_apc$n0[index_1222b_apc], 
                                    log10nmax = param_estimates_apc$nmax[index_1222b_apc], 
                                    mumax = (param_estimates_apc$mumax[index_1222b_apc]*2.303)),
                                  lower = c(0, 0, 0))

growth_parameters_apc$n0[50] <- summary(dat1222b_buchanan_no_lag_apc)$coefficient[1]
growth_parameters_apc$mumax[50] <- summary(dat1222b_buchanan_no_lag_apc)$coefficient[3]
growth_parameters_apc$nmax[50] <- summary(dat1222b_buchanan_no_lag_apc)$coefficient[2]
growth_parameters_apc$aic[50] <- AIC(dat1222b_buchanan_no_lag_apc)
growth_parameters_apc$bic[50] <- BIC(dat1222b_buchanan_no_lag_apc)
growth_parameters_apc$fit[50] <- "yes"
growth_parameters_apc$sampling_code[50] <- "1222_B"
growth_parameters_apc$test[50] <- "apc"
growth_parameters_apc$model[50] <-"buchanan_no_lag"

##----1221A: PC----

index_1221a_pc <- which(param_estimates_pc$sampling_code == "1221_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat1221a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat1221a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_1221a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_1221a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_1221a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[1] <- summary(dat1221a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[1] <- summary(dat1221a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[1] <- summary(dat1221a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[1] <- AIC(dat1221a_baranyi_no_lag_pc)
growth_parameters_pc$bic[1] <- BIC(dat1221a_baranyi_no_lag_pc)
growth_parameters_pc$fit[1] <- "yes"
growth_parameters_pc$sampling_code[1] <- "1221_A"
growth_parameters_pc$test[1] <- "pc"
growth_parameters_pc$model[1] <-"baranyi_no_lag"



#buchanan no lag
dat1221a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat1221a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_1221a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_1221a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_1221a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[2] <- summary(dat1221a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[2] <- summary(dat1221a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[2] <- summary(dat1221a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[2] <- AIC(dat1221a_buchanan_no_lag_pc)
growth_parameters_pc$bic[2] <- BIC(dat1221a_buchanan_no_lag_pc)
growth_parameters_pc$fit[2] <- "yes"
growth_parameters_pc$sampling_code[2] <- "1221_A"
growth_parameters_pc$test[2] <- "pc"
growth_parameters_pc$model[2] <-"buchanan_no_lag"



##----0122A: PC----

index_0122a_pc <- which(param_estimates_pc$sampling_code == "0122_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0122a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0122a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0122a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0122a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0122a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[3] <- summary(dat0122a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[3] <- summary(dat0122a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[3] <- summary(dat0122a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[3] <- AIC(dat0122a_baranyi_no_lag_pc)
growth_parameters_pc$bic[3] <- BIC(dat0122a_baranyi_no_lag_pc)
growth_parameters_pc$fit[3] <- "yes"
growth_parameters_pc$sampling_code[3] <- "0122_A"
growth_parameters_pc$test[3] <- "pc"
growth_parameters_pc$model[3] <-"baranyi_no_lag"



#buchanan no lag
dat0122a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0122a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0122a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0122a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0122a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[4] <- summary(dat0122a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[4] <- summary(dat0122a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[4] <- summary(dat0122a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[4] <- AIC(dat0122a_buchanan_no_lag_pc)
growth_parameters_pc$bic[4] <- BIC(dat0122a_buchanan_no_lag_pc)
growth_parameters_pc$fit[4] <- "yes"
growth_parameters_pc$sampling_code[4] <- "0122_A"
growth_parameters_pc$test[4] <- "pc"
growth_parameters_pc$model[4] <-"buchanan_no_lag"

##----0122B: PC----

index_0122b_pc <- which(param_estimates_pc$sampling_code == "0122_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0122b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0122b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0122b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0122b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0122b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[5] <- summary(dat0122b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[5] <- summary(dat0122b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[5] <- summary(dat0122b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[5] <- AIC(dat0122b_baranyi_no_lag_pc)
growth_parameters_pc$bic[5] <- BIC(dat0122b_baranyi_no_lag_pc)
growth_parameters_pc$fit[5] <- "yes"
growth_parameters_pc$sampling_code[5] <- "0122_B"
growth_parameters_pc$test[5] <- "pc"
growth_parameters_pc$model[5] <-"baranyi_no_lag"



#buchanan no lag
dat0122b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0122b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0122b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0122b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0122b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[6] <- summary(dat0122b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[6] <- summary(dat0122b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[6] <- summary(dat0122b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[6] <- AIC(dat0122b_buchanan_no_lag_pc)
growth_parameters_pc$bic[6] <- BIC(dat0122b_buchanan_no_lag_pc)
growth_parameters_pc$fit[6] <- "yes"
growth_parameters_pc$sampling_code[6] <- "0122_B"
growth_parameters_pc$test[6] <- "pc"
growth_parameters_pc$model[6] <-"buchanan_no_lag"

##----0222A: PC----

index_0222a_pc <- which(param_estimates_pc$sampling_code == "0222_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0222a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0222a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0222a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0222a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0222a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[7] <- summary(dat0222a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[7] <- summary(dat0222a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[7] <- summary(dat0222a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[7] <- AIC(dat0222a_baranyi_no_lag_pc)
growth_parameters_pc$bic[7] <- BIC(dat0222a_baranyi_no_lag_pc)
growth_parameters_pc$fit[7] <- "yes"
growth_parameters_pc$sampling_code[7] <- "0222_A"
growth_parameters_pc$test[7] <- "pc"
growth_parameters_pc$model[7] <-"baranyi_no_lag"



#buchanan no lag
dat0222a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0222a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0222a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0222a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0222a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[8] <- summary(dat0222a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[8] <- summary(dat0222a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[8] <- summary(dat0222a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[8] <- AIC(dat0222a_buchanan_no_lag_pc)
growth_parameters_pc$bic[8] <- BIC(dat0222a_buchanan_no_lag_pc)
growth_parameters_pc$fit[8] <- "yes"
growth_parameters_pc$sampling_code[8] <- "0222_A"
growth_parameters_pc$test[8] <- "pc"
growth_parameters_pc$model[8] <-"buchanan_no_lag"


##----0222B: PC----

index_0222b_pc <- which(param_estimates_pc$sampling_code == "0222_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0222b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0222b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0222b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0222b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0222b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[9] <- summary(dat0222b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[9] <- summary(dat0222b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[9] <- summary(dat0222b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[9] <- AIC(dat0222b_baranyi_no_lag_pc)
growth_parameters_pc$bic[9] <- BIC(dat0222b_baranyi_no_lag_pc)
growth_parameters_pc$fit[9] <- "yes"
growth_parameters_pc$sampling_code[9] <- "0222_B"
growth_parameters_pc$test[9] <- "pc"
growth_parameters_pc$model[9] <-"baranyi_no_lag"



#buchanan no lag
dat0222b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0222b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0222b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0222b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0222b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[10] <- summary(dat0222b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[10] <- summary(dat0222b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[10] <- summary(dat0222b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[10] <- AIC(dat0222b_buchanan_no_lag_pc)
growth_parameters_pc$bic[10] <- BIC(dat0222b_buchanan_no_lag_pc)
growth_parameters_pc$fit[10] <- "yes"
growth_parameters_pc$sampling_code[10] <- "0222_B"
growth_parameters_pc$test[10] <- "pc"
growth_parameters_pc$model[10] <-"buchanan_no_lag"

##----0322A: PC----

index_0322a_pc <- which(param_estimates_pc$sampling_code == "0322_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0322a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0322a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0322a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0322a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0322a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[11] <- summary(dat0322a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[11] <- summary(dat0322a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[11] <- summary(dat0322a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[11] <- AIC(dat0322a_baranyi_no_lag_pc)
growth_parameters_pc$bic[11] <- BIC(dat0322a_baranyi_no_lag_pc)
growth_parameters_pc$fit[11] <- "yes"
growth_parameters_pc$sampling_code[11] <- "0322_A"
growth_parameters_pc$test[11] <- "pc"
growth_parameters_pc$model[11] <-"baranyi_no_lag"



#buchanan no lag
dat0322a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0322a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0322a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0322a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0322a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[12] <- summary(dat0322a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[12] <- summary(dat0322a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[12] <- summary(dat0322a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[12] <- AIC(dat0322a_buchanan_no_lag_pc)
growth_parameters_pc$bic[12] <- BIC(dat0322a_buchanan_no_lag_pc)
growth_parameters_pc$fit[12] <- "yes"
growth_parameters_pc$sampling_code[12] <- "0322_A"
growth_parameters_pc$test[12] <- "pc"
growth_parameters_pc$model[12] <-"buchanan_no_lag"

##----0322B: PC----

index_0322b_pc <- which(param_estimates_pc$sampling_code == "0322_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0322b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0322b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0322b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0322b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0322b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[13] <- summary(dat0322b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[13] <- summary(dat0322b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[13] <- summary(dat0322b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[13] <- AIC(dat0322b_baranyi_no_lag_pc)
growth_parameters_pc$bic[13] <- BIC(dat0322b_baranyi_no_lag_pc)
growth_parameters_pc$fit[13] <- "yes"
growth_parameters_pc$sampling_code[13] <- "0322_B"
growth_parameters_pc$test[13] <- "pc"
growth_parameters_pc$model[13] <-"baranyi_no_lag"



#buchanan no lag
dat0322b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0322b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0322b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0322b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0322b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[14] <- summary(dat0322b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[14] <- summary(dat0322b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[14] <- summary(dat0322b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[14] <- AIC(dat0322b_buchanan_no_lag_pc)
growth_parameters_pc$bic[14] <- BIC(dat0322b_buchanan_no_lag_pc)
growth_parameters_pc$fit[14] <- "yes"
growth_parameters_pc$sampling_code[14] <- "0322_B"
growth_parameters_pc$test[14] <- "pc"
growth_parameters_pc$model[14] <-"buchanan_no_lag"

##----0422A: PC----

index_0422a_pc <- which(param_estimates_pc$sampling_code == "0422_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0422a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0422a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0422a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0422a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0422a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[15] <- summary(dat0422a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[15] <- summary(dat0422a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[15] <- summary(dat0422a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[15] <- AIC(dat0422a_baranyi_no_lag_pc)
growth_parameters_pc$bic[15] <- BIC(dat0422a_baranyi_no_lag_pc)
growth_parameters_pc$fit[15] <- "yes"
growth_parameters_pc$sampling_code[15] <- "0422_A"
growth_parameters_pc$test[15] <- "pc"
growth_parameters_pc$model[15] <-"baranyi_no_lag"



#buchanan no lag
dat0422a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0422a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0422a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0422a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0422a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[16] <- summary(dat0422a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[16] <- summary(dat0422a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[16] <- summary(dat0422a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[16] <- AIC(dat0422a_buchanan_no_lag_pc)
growth_parameters_pc$bic[16] <- BIC(dat0422a_buchanan_no_lag_pc)
growth_parameters_pc$fit[16] <- "yes"
growth_parameters_pc$sampling_code[16] <- "0422_A"
growth_parameters_pc$test[16] <- "pc"
growth_parameters_pc$model[16] <-"buchanan_no_lag"

##----0422B: PC----

index_0422b_pc <- which(param_estimates_pc$sampling_code == "0422_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0422b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0422b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0422b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0422b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0422b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[17] <- summary(dat0422b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[17] <- summary(dat0422b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[17] <- summary(dat0422b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[17] <- AIC(dat0422b_baranyi_no_lag_pc)
growth_parameters_pc$bic[17] <- BIC(dat0422b_baranyi_no_lag_pc)
growth_parameters_pc$fit[17] <- "no"
growth_parameters_pc$sampling_code[17] <- "0422_B"
growth_parameters_pc$test[17] <- "pc"
growth_parameters_pc$model[17] <-"baranyi_no_lag"



#buchanan no lag
dat0422b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0422b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0422b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0422b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0422b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[18] <- summary(dat0422b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[18] <- summary(dat0422b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[18] <- summary(dat0422b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[18] <- AIC(dat0422b_buchanan_no_lag_pc)
growth_parameters_pc$bic[18] <- BIC(dat0422b_buchanan_no_lag_pc)
growth_parameters_pc$fit[18] <- "no"
growth_parameters_pc$sampling_code[18] <- "0422_B"
growth_parameters_pc$test[18] <- "pc"
growth_parameters_pc$model[18] <-"buchanan_no_lag"

##----0522A: PC----

index_0522a_pc <- which(param_estimates_pc$sampling_code == "0522_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0522a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0522a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0522a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0522a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0522a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[19] <- summary(dat0522a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[19] <- summary(dat0522a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[19] <- summary(dat0522a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[19] <- AIC(dat0522a_baranyi_no_lag_pc)
growth_parameters_pc$bic[19] <- BIC(dat0522a_baranyi_no_lag_pc)
growth_parameters_pc$fit[19] <- "no"
growth_parameters_pc$sampling_code[19] <- "0522_A"
growth_parameters_pc$test[19] <- "pc"
growth_parameters_pc$model[19] <-"baranyi_no_lag"



#buchanan no lag
dat0522a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0522a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0522a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0522a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0522a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[20] <- summary(dat0522a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[20] <- summary(dat0522a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[20] <- summary(dat0522a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[20] <- AIC(dat0522a_buchanan_no_lag_pc)
growth_parameters_pc$bic[20] <- BIC(dat0522a_buchanan_no_lag_pc)
growth_parameters_pc$fit[20] <- "no"
growth_parameters_pc$sampling_code[20] <- "0522_A"
growth_parameters_pc$test[20] <- "pc"
growth_parameters_pc$model[20] <-"buchanan_no_lag"

##----0522B: PC----

index_0522b_pc <- which(param_estimates_pc$sampling_code == "0522_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0522b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0522b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0522b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0522b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0522b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[21] <- summary(dat0522b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[21] <- summary(dat0522b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[21] <- summary(dat0522b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[21] <- AIC(dat0522b_baranyi_no_lag_pc)
growth_parameters_pc$bic[21] <- BIC(dat0522b_baranyi_no_lag_pc)
growth_parameters_pc$fit[21] <- "yes"
growth_parameters_pc$sampling_code[21] <- "0522_B"
growth_parameters_pc$test[21] <- "pc"
growth_parameters_pc$model[21] <-"baranyi_no_lag"



#buchanan no lag
dat0522b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0522b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0522b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0522b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0522b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[22] <- summary(dat0522b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[22] <- summary(dat0522b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[22] <- summary(dat0522b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[22] <- AIC(dat0522b_buchanan_no_lag_pc)
growth_parameters_pc$bic[22] <- BIC(dat0522b_buchanan_no_lag_pc)
growth_parameters_pc$fit[22] <- "yes"
growth_parameters_pc$sampling_code[22] <- "0522_B"
growth_parameters_pc$test[22] <- "pc"
growth_parameters_pc$model[22] <-"buchanan_no_lag"

##----0522C: PC----

index_0522c_pc <- which(param_estimates_pc$sampling_code == "0522_C" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0522c_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0522c_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0522c_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0522c_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0522c_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[23] <- summary(dat0522c_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[23] <- summary(dat0522c_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[23] <- summary(dat0522c_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[23] <- AIC(dat0522c_baranyi_no_lag_pc)
growth_parameters_pc$bic[23] <- BIC(dat0522c_baranyi_no_lag_pc)
growth_parameters_pc$fit[23] <- "yes"
growth_parameters_pc$sampling_code[23] <- "0522_C"
growth_parameters_pc$test[23] <- "pc"
growth_parameters_pc$model[23] <-"baranyi_no_lag"



#buchanan no lag
dat0522c_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0522c_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0522c_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0522c_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0522c_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[24] <- summary(dat0522c_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[24] <- summary(dat0522c_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[24] <- summary(dat0522c_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[24] <- AIC(dat0522c_buchanan_no_lag_pc)
growth_parameters_pc$bic[24] <- BIC(dat0522c_buchanan_no_lag_pc)
growth_parameters_pc$fit[24] <- "yes"
growth_parameters_pc$sampling_code[24] <- "0522_C"
growth_parameters_pc$test[24] <- "pc"
growth_parameters_pc$model[24] <-"buchanan_no_lag"

##----0622A: PC----

index_0622a_pc <- which(param_estimates_pc$sampling_code == "0622_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0622a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0622a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0622a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0622a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0622a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[25] <- summary(dat0622a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[25] <- summary(dat0622a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[25] <- summary(dat0622a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[25] <- AIC(dat0622a_baranyi_no_lag_pc)
growth_parameters_pc$bic[25] <- BIC(dat0622a_baranyi_no_lag_pc)
growth_parameters_pc$fit[25] <- "yes"
growth_parameters_pc$sampling_code[25] <- "0622_A"
growth_parameters_pc$test[25] <- "pc"
growth_parameters_pc$model[25] <-"baranyi_no_lag"



#buchanan no lag
dat0622a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0622a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0622a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0622a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0622a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[26] <- summary(dat0622a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[26] <- summary(dat0622a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[26] <- summary(dat0622a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[26] <- AIC(dat0622a_buchanan_no_lag_pc)
growth_parameters_pc$bic[26] <- BIC(dat0622a_buchanan_no_lag_pc)
growth_parameters_pc$fit[26] <- "yes"
growth_parameters_pc$sampling_code[26] <- "0622_A"
growth_parameters_pc$test[26] <- "pc"
growth_parameters_pc$model[26] <-"buchanan_no_lag"

##----0622B: PC----

index_0622b_pc <- which(param_estimates_pc$sampling_code == "0622_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0622b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0622b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0622b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0622b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0622b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[27] <- summary(dat0622b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[27] <- summary(dat0622b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[27] <- summary(dat0622b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[27] <- AIC(dat0622b_baranyi_no_lag_pc)
growth_parameters_pc$bic[27] <- BIC(dat0622b_baranyi_no_lag_pc)
growth_parameters_pc$fit[27] <- "yes"
growth_parameters_pc$sampling_code[27] <- "0622_B"
growth_parameters_pc$test[27] <- "pc"
growth_parameters_pc$model[27] <-"baranyi_no_lag"



#buchanan no lag
dat0622b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0622b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0622b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0622b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0622b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[28] <- summary(dat0622b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[28] <- summary(dat0622b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[28] <- summary(dat0622b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[28] <- AIC(dat0622b_buchanan_no_lag_pc)
growth_parameters_pc$bic[28] <- BIC(dat0622b_buchanan_no_lag_pc)
growth_parameters_pc$fit[28] <- "yes"
growth_parameters_pc$sampling_code[28] <- "0622_B"
growth_parameters_pc$test[28] <- "pc"
growth_parameters_pc$model[28] <-"buchanan_no_lag"

##----0722A: PC----

index_0722a_pc <- which(param_estimates_pc$sampling_code == "0722_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0722a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0722a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0722a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0722a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0722a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[29] <- summary(dat0722a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[29] <- summary(dat0722a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[29] <- summary(dat0722a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[29] <- AIC(dat0722a_baranyi_no_lag_pc)
growth_parameters_pc$bic[29] <- BIC(dat0722a_baranyi_no_lag_pc)
growth_parameters_pc$fit[29] <- "yes"
growth_parameters_pc$sampling_code[29] <- "0722_A"
growth_parameters_pc$test[29] <- "pc"
growth_parameters_pc$model[29] <-"baranyi_no_lag"



#buchanan no lag
dat0722a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0722a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0722a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0722a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0722a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[30] <- summary(dat0722a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[30] <- summary(dat0722a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[30] <- summary(dat0722a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[30] <- AIC(dat0722a_buchanan_no_lag_pc)
growth_parameters_pc$bic[30] <- BIC(dat0722a_buchanan_no_lag_pc)
growth_parameters_pc$fit[30] <- "yes"
growth_parameters_pc$sampling_code[30] <- "0722_A"
growth_parameters_pc$test[30] <- "pc"
growth_parameters_pc$model[30] <-"buchanan_no_lag"

##----0822A: PC----

index_0822a_pc <- which(param_estimates_pc$sampling_code == "0822_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0822a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0822a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0822a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0822a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0822a_pc]*2.303)),
                                     lower = c(0, 0, 0),
                                     control = nls.lm.control(maxiter = 100))

growth_parameters_pc$n0[31] <- summary(dat0822a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[31] <- summary(dat0822a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[31] <- summary(dat0822a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[31] <- AIC(dat0822a_baranyi_no_lag_pc)
growth_parameters_pc$bic[31] <- BIC(dat0822a_baranyi_no_lag_pc)
growth_parameters_pc$fit[31] <- "yes"
growth_parameters_pc$sampling_code[31] <- "0822_A"
growth_parameters_pc$test[31] <- "pc"
growth_parameters_pc$model[31] <-"baranyi_no_lag"



#buchanan no lag
dat0822a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0822a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0822a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0822a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0822a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[32] <- summary(dat0822a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[32] <- summary(dat0822a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[32] <- summary(dat0822a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[32] <- AIC(dat0822a_buchanan_no_lag_pc)
growth_parameters_pc$bic[32] <- BIC(dat0822a_buchanan_no_lag_pc)
growth_parameters_pc$fit[32] <- "yes"
growth_parameters_pc$sampling_code[32] <- "0822_A"
growth_parameters_pc$test[32] <- "pc"
growth_parameters_pc$model[32] <-"buchanan_no_lag"

##----0822B: PC----

index_0822b_pc <- which(param_estimates_pc$sampling_code == "0822_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0822b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0822b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0822b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0822b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0822b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[33] <- summary(dat0822b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[33] <- summary(dat0822b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[33] <- summary(dat0822b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[33] <- AIC(dat0822b_baranyi_no_lag_pc)
growth_parameters_pc$bic[33] <- BIC(dat0822b_baranyi_no_lag_pc)
growth_parameters_pc$fit[33] <- "yes"
growth_parameters_pc$sampling_code[33] <- "0822_B"
growth_parameters_pc$test[33] <- "pc"
growth_parameters_pc$model[33] <-"baranyi_no_lag"



#buchanan no lag
dat0822b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0822b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0822b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0822b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0822b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[34] <- summary(dat0822b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[34] <- summary(dat0822b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[34] <- summary(dat0822b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[34] <- AIC(dat0822b_buchanan_no_lag_pc)
growth_parameters_pc$bic[34] <- BIC(dat0822b_buchanan_no_lag_pc)
growth_parameters_pc$fit[34] <- "yes"
growth_parameters_pc$sampling_code[34] <- "0822_B"
growth_parameters_pc$test[34] <- "pc"
growth_parameters_pc$model[34] <-"buchanan_no_lag"

##----0922A: PC----

index_0922a_pc <- which(param_estimates_pc$sampling_code == "0922_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0922a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0922a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0922a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0922a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0922a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[35] <- summary(dat0922a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[35] <- summary(dat0922a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[35] <- summary(dat0922a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[35] <- AIC(dat0922a_baranyi_no_lag_pc)
growth_parameters_pc$bic[35] <- BIC(dat0922a_baranyi_no_lag_pc)
growth_parameters_pc$fit[35] <- "yes"
growth_parameters_pc$sampling_code[35] <- "0922_A"
growth_parameters_pc$test[35] <- "pc"
growth_parameters_pc$model[35] <-"baranyi_no_lag"



#buchanan no lag
dat0922a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0922a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0922a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0922a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0922a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[36] <- summary(dat0922a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[36] <- summary(dat0922a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[36] <- summary(dat0922a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[36] <- AIC(dat0922a_buchanan_no_lag_pc)
growth_parameters_pc$bic[36] <- BIC(dat0922a_buchanan_no_lag_pc)
growth_parameters_pc$fit[36] <- "yes"
growth_parameters_pc$sampling_code[36] <- "0922_A"
growth_parameters_pc$test[36] <- "pc"
growth_parameters_pc$model[36] <-"buchanan_no_lag"

##----0922B: PC----

index_0922b_pc <- which(param_estimates_pc$sampling_code == "0922_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat0922b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat0922b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_0922b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_0922b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_0922b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[37] <- summary(dat0922b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[37] <- summary(dat0922b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[37] <- summary(dat0922b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[37] <- AIC(dat0922b_baranyi_no_lag_pc)
growth_parameters_pc$bic[37] <- BIC(dat0922b_baranyi_no_lag_pc)
growth_parameters_pc$fit[37] <- "yes"
growth_parameters_pc$sampling_code[37] <- "0922_B"
growth_parameters_pc$test[37] <- "pc"
growth_parameters_pc$model[37] <-"baranyi_no_lag"

#buchanan no lag
dat0922b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat0922b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_0922b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_0922b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_0922b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[38] <- summary(dat0922b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[38] <- summary(dat0922b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[38] <- summary(dat0922b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[38] <- AIC(dat0922b_buchanan_no_lag_pc)
growth_parameters_pc$bic[38] <- BIC(dat0922b_buchanan_no_lag_pc)
growth_parameters_pc$fit[38] <- "yes"
growth_parameters_pc$sampling_code[38] <- "0922_B"
growth_parameters_pc$test[38] <- "pc"
growth_parameters_pc$model[38] <-"buchanan_no_lag"

##----1022A: PC----

index_1022a_pc <- which(param_estimates_pc$sampling_code == "1022_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat1022a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat1022a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_1022a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_1022a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_1022a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[39] <- summary(dat1022a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[39] <- summary(dat1022a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[39] <- summary(dat1022a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[39] <- AIC(dat1022a_baranyi_no_lag_pc)
growth_parameters_pc$bic[39] <- BIC(dat1022a_baranyi_no_lag_pc)
growth_parameters_pc$fit[39] <- "yes"
growth_parameters_pc$sampling_code[39] <- "1022_A"
growth_parameters_pc$test[39] <- "pc"
growth_parameters_pc$model[39] <-"baranyi_no_lag"

#buchanan no lag
dat1022a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat1022a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_1022a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_1022a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_1022a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[40] <- summary(dat1022a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[40] <- summary(dat1022a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[40] <- summary(dat1022a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[40] <- AIC(dat1022a_buchanan_no_lag_pc)
growth_parameters_pc$bic[40] <- BIC(dat1022a_buchanan_no_lag_pc)
growth_parameters_pc$fit[40] <- "yes"
growth_parameters_pc$sampling_code[40] <- "1022_A"
growth_parameters_pc$test[40] <- "pc"
growth_parameters_pc$model[40] <-"buchanan_no_lag"

##----1022B: PC----

index_1022b_pc <- which(param_estimates_pc$sampling_code == "1022_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat1022b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat1022b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_1022b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_1022b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_1022b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[41] <- summary(dat1022b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[41] <- summary(dat1022b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[41] <- summary(dat1022b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[41] <- AIC(dat1022b_baranyi_no_lag_pc)
growth_parameters_pc$bic[41] <- BIC(dat1022b_baranyi_no_lag_pc)
growth_parameters_pc$fit[41] <- "yes"
growth_parameters_pc$sampling_code[41] <- "1022_B"
growth_parameters_pc$test[41] <- "pc"
growth_parameters_pc$model[41] <-"baranyi_no_lag"

#buchanan no lag
dat1022b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat1022b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_1022b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_1022b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_1022b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[42] <- summary(dat1022b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[42] <- summary(dat1022b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[42] <- summary(dat1022b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[42] <- AIC(dat1022b_buchanan_no_lag_pc)
growth_parameters_pc$bic[42] <- BIC(dat1022b_buchanan_no_lag_pc)
growth_parameters_pc$fit[42] <- "yes"
growth_parameters_pc$sampling_code[42] <- "1022_B"
growth_parameters_pc$test[42] <- "pc"
growth_parameters_pc$model[42] <-"buchanan_no_lag"

##----1122A: PC----

index_1122a_pc <- which(param_estimates_pc$sampling_code == "1122_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat1122a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat1122a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_1122a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_1122a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_1122a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[43] <- summary(dat1122a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[43] <- summary(dat1122a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[43] <- summary(dat1122a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[43] <- AIC(dat1122a_baranyi_no_lag_pc)
growth_parameters_pc$bic[43] <- BIC(dat1122a_baranyi_no_lag_pc)
growth_parameters_pc$fit[43] <- "yes"
growth_parameters_pc$sampling_code[43] <- "1122_A"
growth_parameters_pc$test[43] <- "pc"
growth_parameters_pc$model[43] <-"baranyi_no_lag"

#buchanan no lag
dat1122a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat1122a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_1122a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_1122a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_1122a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[44] <- summary(dat1122a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[44] <- summary(dat1122a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[44] <- summary(dat1122a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[44] <- AIC(dat1122a_buchanan_no_lag_pc)
growth_parameters_pc$bic[44] <- BIC(dat1122a_buchanan_no_lag_pc)
growth_parameters_pc$fit[44] <- "yes"
growth_parameters_pc$sampling_code[44] <- "1122_A"
growth_parameters_pc$test[44] <- "pc"
growth_parameters_pc$model[44] <-"buchanan_no_lag"

##----1122B: PC----

index_1122b_pc <- which(param_estimates_pc$sampling_code == "1122_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat1122b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat1122b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_1122b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_1122b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_1122b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[45] <- summary(dat1122b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[45] <- summary(dat1122b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[45] <- summary(dat1122b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[45] <- AIC(dat1122b_baranyi_no_lag_pc)
growth_parameters_pc$bic[45] <- BIC(dat1122b_baranyi_no_lag_pc)
growth_parameters_pc$fit[45] <- "yes"
growth_parameters_pc$sampling_code[45] <- "1122_B"
growth_parameters_pc$test[45] <- "pc"
growth_parameters_pc$model[45] <-"baranyi_no_lag"

#buchanan no lag
dat1122b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat1122b_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_1122b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_1122b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_1122b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[46] <- summary(dat1122b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[46] <- summary(dat1122b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[46] <- summary(dat1122b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[46] <- AIC(dat1122b_buchanan_no_lag_pc)
growth_parameters_pc$bic[46] <- BIC(dat1122b_buchanan_no_lag_pc)
growth_parameters_pc$fit[46] <- "yes"
growth_parameters_pc$sampling_code[46] <- "1122_B"
growth_parameters_pc$test[46] <- "pc"
growth_parameters_pc$model[46] <-"buchanan_no_lag"

##----1222A: PC----

index_1222a_pc <- which(param_estimates_pc$sampling_code == "1222_A" & param_estimates_pc$test == "pc")

#baranyi no lag
dat1222a_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat1222a_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_1222a_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_1222a_pc], 
                                       mumax = (param_estimates_pc$mumax[index_1222a_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[47] <- summary(dat1222a_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[47] <- summary(dat1222a_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[47] <- summary(dat1222a_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[47] <- AIC(dat1222a_baranyi_no_lag_pc)
growth_parameters_pc$bic[47] <- BIC(dat1222a_baranyi_no_lag_pc)
growth_parameters_pc$fit[47] <- "yes"
growth_parameters_pc$sampling_code[47] <- "1222_A"
growth_parameters_pc$test[47] <- "pc"
growth_parameters_pc$model[47] <-"baranyi_no_lag"

#buchanan no lag
dat1222a_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat1222a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_1222a_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_1222a_pc], 
                                        mumax = (param_estimates_pc$mumax[index_1222a_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[48] <- summary(dat1222a_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[48] <- summary(dat1222a_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[48] <- summary(dat1222a_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[48] <- AIC(dat1222a_buchanan_no_lag_pc)
growth_parameters_pc$bic[48] <- BIC(dat1222a_buchanan_no_lag_pc)
growth_parameters_pc$fit[48] <- "yes"
growth_parameters_pc$sampling_code[48] <- "1222_A"
growth_parameters_pc$test[48] <- "pc"
growth_parameters_pc$model[48] <-"buchanan_no_lag"

##----1222B: PC----

index_1222b_pc <- which(param_estimates_pc$sampling_code == "1222_B" & param_estimates_pc$test == "pc")

#baranyi no lag
dat1222b_baranyi_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                       baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = dat1222b_pc, 
                                     start=list(
                                       log10n0 = param_estimates_pc$n0[index_1222b_pc], 
                                       log10nmax = param_estimates_pc$nmax[index_1222b_pc], 
                                       mumax = (param_estimates_pc$mumax[index_1222b_pc]*2.303)),
                                     lower = c(0, 0, 0))

growth_parameters_pc$n0[49] <- summary(dat1222b_baranyi_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[49] <- summary(dat1222b_baranyi_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[49] <- summary(dat1222b_baranyi_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[49] <- AIC(dat1222b_baranyi_no_lag_pc)
growth_parameters_pc$bic[49] <- BIC(dat1222b_baranyi_no_lag_pc)
growth_parameters_pc$fit[49] <- "yes"
growth_parameters_pc$sampling_code[49] <- "1222_B"
growth_parameters_pc$test[49] <- "pc"
growth_parameters_pc$model[49] <-"baranyi_no_lag"

#buchanan no lag
dat1222b_buchanan_no_lag_pc <- nlsLM(log_conc_geom ~ 
                                        buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = dat1222a_pc, 
                                      start=list(
                                        log10n0 = param_estimates_pc$n0[index_1222b_pc], 
                                        log10nmax = param_estimates_pc$nmax[index_1222b_pc], 
                                        mumax = (param_estimates_pc$mumax[index_1222b_pc]*2.303)),
                                      lower = c(0, 0, 0))

growth_parameters_pc$n0[50] <- summary(dat1222b_buchanan_no_lag_pc)$coefficient[1]
growth_parameters_pc$mumax[50] <- summary(dat1222b_buchanan_no_lag_pc)$coefficient[3]
growth_parameters_pc$nmax[50] <- summary(dat1222b_buchanan_no_lag_pc)$coefficient[2]
growth_parameters_pc$aic[50] <- AIC(dat1222b_buchanan_no_lag_pc)
growth_parameters_pc$bic[50] <- BIC(dat1222b_buchanan_no_lag_pc)
growth_parameters_pc$fit[50] <- "yes"
growth_parameters_pc$sampling_code[50] <- "1222_B"
growth_parameters_pc$test[50] <- "pc"
growth_parameters_pc$model[50] <-"buchanan_no_lag"

## -----------------------------Selecting the model with the best fit-----------

best_fit_apc <- growth_parameters_apc %>%
  filter(fit == "yes") %>%
  group_by(sampling_code) %>%
  arrange(aic) %>%
  slice(1)

best_fit_pc <- growth_parameters_pc %>%
  filter(fit == "yes") %>%
  group_by(sampling_code) %>%
  arrange(aic) %>%
  slice(1)

best_fit <- bind_rows(best_fit_apc, best_fit_pc)

## -----------------------------Exporting data----------------------------------

write.csv(best_fit, "output/primary_growth_model_parameters.csv", row.names = FALSE)

