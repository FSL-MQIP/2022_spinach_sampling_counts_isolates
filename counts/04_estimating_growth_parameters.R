## ------------------------------Title------------------------------------------
# Estimating initial parameters in order to fit growth curves to wrangled data from the 2022 sampling (data collected up till 3 January 2022)

## ------------------------------Description------------------------------------
# Project: CIDA Spinach 

# Script description: Estimating growth parameters for fitting models without lag (Baranyi no lag and Buchanan no lag). 
#                     Models without lag will be used for fitting, since it appears that there is no apparent lag phase (based on plots of microbial concentration over shelf life)
#                     Growth parameters will be estimated based on changes in APC of the samples

## ------------------------------Packages---------------------------------------
# Loading packages
library(tidyverse); library(dplyr);

## ------------------------------Raw Data---------------------------------------

conc <- read.csv("data/wrangled/shelf_life_conc_by_testing_day.csv", header = TRUE)

## ------------------------------Wrangling--------------------------------------
#Filtering data by sampling and shelf life

dat1221a_apc <- filter(conc, sampling_code == "1221_A" & test == "APC" & day != "H")

dat0122a_apc <- filter(conc, sampling_code == "0122_A" & test == "APC" & day != "H")
dat0122b_apc <- filter(conc, sampling_code == "0122_B" & test == "APC" & day != "H")

dat0222a_apc <- filter(conc, sampling_code == "0222_A" & test == "APC" & day != "H")
dat0222b_apc <- filter(conc, sampling_code == "0222_B" & test == "APC" & day != "H")

dat0322a_apc <- filter(conc, sampling_code == "0322_A" & test == "APC" & day != "H")
dat0322b_apc <- filter(conc, sampling_code == "0322_B" & test == "APC" & day != "H")

dat0422a_apc <- filter(conc, sampling_code == "0422_A" & test == "APC" & day != "H")
dat0422b_apc <- filter(conc, sampling_code == "0422_B" & test == "APC" & day != "H")

dat0522a_apc <- filter(conc, sampling_code == "0522_A" & test == "APC" & day != "H")
dat0522b_apc <- filter(conc, sampling_code == "0522_B" & test == "APC" & day != "H")
dat0522c_apc <- filter(conc, sampling_code == "0522_C" & test == "APC" & day != "H")

dat0622a_apc <- filter(conc, sampling_code == "0622_A" & test == "APC" & day != "H")
dat0622b_apc <- filter(conc, sampling_code == "0622_B" & test == "APC" & day != "H")

dat0722a_apc <- filter(conc, sampling_code == "0722_A" & test == "APC" & day != "H")

dat0822a_apc <- filter(conc, sampling_code == "0822_A" & test == "APC" & day != "H")
dat0822b_apc <- filter(conc, sampling_code == "0822_B" & test == "APC" & day != "H")

dat0922a_apc <- filter(conc, sampling_code == "0922_A" & test == "APC" & day != "H")
dat0922b_apc <- filter(conc, sampling_code == "0922_B" & test == "APC" & day != "H")

dat1022a_apc <- filter(conc, sampling_code == "1022_A" & test == "APC" & day != "H")
dat1022b_apc <- filter(conc, sampling_code == "1022_B" & test == "APC" & day != "H")

dat1122a_apc <- filter(conc, sampling_code == "1122_A" & test == "APC" & day != "H")
dat1122b_apc <- filter(conc, sampling_code == "1122_B" & test == "APC" & day != "H")

dat1222a_apc <- filter(conc, sampling_code == "1222_A" & test == "APC" & day != "H")
dat1222b_apc <- filter(conc, sampling_code == "1222_B" & test == "APC" & day != "H")

dat1221a_pc <- filter(conc, sampling_code == "1221_A" & test == "PC" & day != "H")

dat0122a_pc <- filter(conc, sampling_code == "0122_A" & test == "PC" & day != "H")
dat0122b_pc <- filter(conc, sampling_code == "0122_B" & test == "PC" & day != "H")

dat0222a_pc <- filter(conc, sampling_code == "0222_A" & test == "PC" & day != "H")
dat0222b_pc <- filter(conc, sampling_code == "0222_B" & test == "PC" & day != "H")

dat0322a_pc <- filter(conc, sampling_code == "0322_A" & test == "PC" & day != "H")
dat0322b_pc <- filter(conc, sampling_code == "0322_B" & test == "PC" & day != "H")

dat0422a_pc <- filter(conc, sampling_code == "0422_A" & test == "PC" & day != "H")
dat0422b_pc <- filter(conc, sampling_code == "0422_B" & test == "PC" & day != "H")

dat0522a_pc <- filter(conc, sampling_code == "0522_A" & test == "PC" & day != "H")
dat0522b_pc <- filter(conc, sampling_code == "0522_B" & test == "PC" & day != "H")
dat0522c_pc <- filter(conc, sampling_code == "0522_C" & test == "PC" & day != "H")

dat0622a_pc <- filter(conc, sampling_code == "0622_A" & test == "PC" & day != "H")
dat0622b_pc <- filter(conc, sampling_code == "0622_B" & test == "PC" & day != "H")

dat0722a_pc <- filter(conc, sampling_code == "0722_A" & test == "PC" & day != "H")

dat0822a_pc <- filter(conc, sampling_code == "0822_A" & test == "PC" & day != "H")
dat0822b_pc <- filter(conc, sampling_code == "0822_B" & test == "PC" & day != "H")

dat0922a_pc <- filter(conc, sampling_code == "0922_A" & test == "PC" & day != "H")
dat0922b_pc <- filter(conc, sampling_code == "0922_B" & test == "PC" & day != "H")

dat1022a_pc <- filter(conc, sampling_code == "1022_A" & test == "PC" & day != "H")
dat1022b_pc <- filter(conc, sampling_code == "1022_B" & test == "PC" & day != "H")

dat1122a_pc <- filter(conc, sampling_code == "1122_A" & test == "PC" & day != "H")
dat1122b_pc <- filter(conc, sampling_code == "1122_B" & test == "PC" & day != "H")

dat1222a_pc <- filter(conc, sampling_code == "1222_A" & test == "PC" & day != "H")
dat1222b_pc <- filter(conc, sampling_code == "1222_B" & test == "PC" & day != "H")

## -----------------------------Estimating growth parameters--------------------
#Setting up a dataframe for recording growth parameters 
initial_parameters_apc <- data.frame(matrix(data = NA, nrow = 25, ncol = 5))
colnames(initial_parameters_apc) <- c("sampling_code", "n0", "mumax", "nmax", "test")

initial_parameters_pc <- data.frame(matrix(data = NA, nrow = 25, ncol = 5))
colnames(initial_parameters_pc) <- c("sampling_code", "n0", "mumax", "nmax", "test")

##----1221A: APC----

initial_parameters_apc[1, 1] <- "1221_A"
initial_parameters_apc[1, 5] <- "apc"

#n0 
initial_parameters_apc[1, 2] <- dat1221a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[1,4] <- dat1221a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1221a_apc <- dat1221a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1221a_apc <- dat1221a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_apc[1,3] <- (initial_parameters_apc[1,4] - initial_parameters_apc[1, 2])/(max_day_1221a_apc - min_day_1221a_apc)

##----1221A: PC----

initial_parameters_pc[1, 1] <- "1221_A"
initial_parameters_pc[1, 5] <- "pc"

#n0 
initial_parameters_pc[1, 2] <- dat1221a_pc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[1,4] <- dat1221a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1221a_pc <- dat1221a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1221a_pc <- dat1221a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_pc[1,3] <- (initial_parameters_pc[1,4] - initial_parameters_pc[1, 2])/(max_day_1221a_pc - min_day_1221a_pc)

##----0122A: APC----

initial_parameters_apc[2, 1] <- "0122_A"
initial_parameters_apc[2, 5] <- "apc"

#n0 
initial_parameters_apc[2, 2] <- dat0122a_apc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[2,4] <- dat0122a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0122a_apc <- dat0122a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0122a_apc <- dat0122a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[2,3] <- (initial_parameters_apc[2,4] - initial_parameters_apc[2, 2])/(max_day_0122a_apc - min_day_0122a_apc)

##----0122A: PC----

initial_parameters_pc[2, 1] <- "0122_A"
initial_parameters_pc[2, 5] <- "pc"

#n0 
initial_parameters_pc[2, 2] <- dat0122a_pc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[2,4] <- dat0122a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0122a_pc <- dat0122a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0122a_pc <- dat0122a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_pc[2,3] <- (initial_parameters_pc[2,4] - initial_parameters_pc[2, 2])/(max_day_0122a_pc - min_day_0122a_pc)

##----0122B: APC----

initial_parameters_apc[3, 1] <- "0122_B"
initial_parameters_apc[3, 5] <- "apc"

#n0 
initial_parameters_apc[3, 2] <- dat0122b_apc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[3,4] <- dat0122b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0122b_apc <- dat0122b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0122b_apc <- dat0122b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_apc[3,3] <- (initial_parameters_apc[3,4] - initial_parameters_apc[3, 2])/(max_day_0122b_apc - min_day_0122b_apc)

##----0122B: PC----

initial_parameters_pc[3, 1] <- "0122_B"
initial_parameters_pc[3, 5] <- "pc"

#n0 
initial_parameters_pc[3, 2] <- dat0122b_pc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[3,4] <- dat0122b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0122b_pc <- dat0122b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0122b_pc <- dat0122b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_pc[3,3] <- (initial_parameters_pc[3,4] - initial_parameters_pc[3, 2])/(max_day_0122b_pc - min_day_0122b_pc)

##----0222A: APC----

initial_parameters_apc[4, 1] <- "0222_A"
initial_parameters_apc[4, 5] <- "apc"

#n0 
initial_parameters_apc[4, 2] <- dat0222a_apc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[4,4] <- dat0222a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0222a_apc <- dat0222a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0222a_apc <- dat0222a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_apc[4,3] <- (initial_parameters_apc[4,4] - initial_parameters_apc[4, 2])/(max_day_0222a_apc - min_day_0222a_apc)

##----0222A: PC----

initial_parameters_pc[4, 1] <- "0222_A"
initial_parameters_pc[4, 5] <- "pc"

#n0 
initial_parameters_pc[4, 2] <- dat0222a_pc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[4,4] <- dat0222a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0222a_pc <- dat0222a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0222a_pc <- dat0222a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_pc[4,3] <- (initial_parameters_pc[4,4] - initial_parameters_pc[4, 2])/(max_day_0222a_pc - min_day_0222a_pc)

##----0222B: APC----

initial_parameters_apc[5, 1] <- "0222_B"
initial_parameters_apc[5, 5] <- "apc"

#n0 
initial_parameters_apc[5, 2] <- dat0222b_apc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[5,4] <- dat0222b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0222b_apc <- dat0222b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0222b_apc <- dat0222b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_apc[5,3] <- (initial_parameters_apc[5,4] - initial_parameters_apc[5, 2])/(max_day_0222b_apc - min_day_0222b_apc)

##----0222B: PC----

initial_parameters_pc[5, 1] <- "0222_B"
initial_parameters_pc[5, 5] <- "pc"

#n0 
initial_parameters_pc[5, 2] <- dat0222b_pc  %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[5,4] <- dat0222b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0222b_pc <- dat0222b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0222b_pc <- dat0222b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_pc[5,3] <- (initial_parameters_pc[5,4] - initial_parameters_pc[5, 2])/(max_day_0222b_pc - min_day_0222b_pc)

##----0322A: APC----

initial_parameters_apc[6, 1] <- "0322_A"
initial_parameters_apc[6, 5] <- "apc"

#n0 
initial_parameters_apc[6, 2] <- dat0322a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[6,4] <- dat0322a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0322a_apc <- dat0322a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0322a_apc <- dat0322a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_apc[6,3] <- (initial_parameters_apc[6,4] - initial_parameters_apc[6, 2])/(max_day_0322a_apc - min_day_0322a_apc)

##----0322A: PC----

initial_parameters_pc[6, 1] <- "0322_A"
initial_parameters_pc[6, 5] <- "pc"

#n0 
initial_parameters_pc[6, 2] <- dat0322a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[6,4] <- dat0322a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0322a_pc <- dat0322a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0322a_pc <- dat0322a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[6,3] <- (initial_parameters_pc[6,4] - initial_parameters_pc[6, 2])/(max_day_0322a_pc - min_day_0322a_pc)

##----0322B: APC----

initial_parameters_apc[7, 1] <- "0322_B"
initial_parameters_apc[7, 5] <- "apc"

#n0 
initial_parameters_apc[7, 2] <- dat0322b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[7,4] <- dat0322b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0322b_apc <- dat0322b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0322b_apc <- dat0322b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[7,3] <- (initial_parameters_apc[7,4] - initial_parameters_apc[7,2])/(max_day_0322b_apc - min_day_0322b_apc)

##----0322B: PC----

initial_parameters_pc[7, 1] <- "0322_B"
initial_parameters_pc[7, 5] <- "pc"

#n0 
initial_parameters_pc[7, 2] <- dat0322b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[7,4] <- dat0322b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0322b_pc <- dat0322b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0322b_pc <- dat0322b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[7,3] <- (initial_parameters_pc[7,4] - initial_parameters_pc[7,2])/(max_day_0322b_pc - min_day_0322b_pc)

##----0422A: APC----

initial_parameters_apc[8, 1] <- "0422_A"
initial_parameters_apc[8, 5] <- "apc"

#n0 
initial_parameters_apc[8, 2] <- dat0422a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[8,4] <- dat0422a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0422a_apc <- dat0422a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0422a_apc <- dat0422a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[8,3] <- (initial_parameters_apc[8,4] - initial_parameters_apc[8,2])/(max_day_0422a_apc - min_day_0422a_apc)

##----0422A: PC----

initial_parameters_pc[8, 1] <- "0422_A"
initial_parameters_pc[8, 5] <- "pc"

#n0 
initial_parameters_pc[8, 2] <- dat0422a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[8,4] <- dat0422a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0422a_pc <- dat0422a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0422a_pc <- dat0422a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[8,3] <- (initial_parameters_pc[8,4] - initial_parameters_pc[8,2])/(max_day_0422a_pc - min_day_0422a_pc)

##----0422B: APC----

initial_parameters_apc[9, 1] <- "0422_B"
initial_parameters_apc[9, 5] <- "apc"

#n0 
initial_parameters_apc[9, 2] <- dat0422b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[9,4] <- dat0422b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0422b_apc <- dat0422b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0422b_apc <- dat0422b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[9,3] <- (initial_parameters_apc[9,4] - initial_parameters_apc[9,2])/(max_day_0422b_apc - min_day_0422b_apc)

##----0422B: PC----

initial_parameters_pc[9, 1] <- "0422_B"
initial_parameters_pc[9, 5] <- "pc"

#n0 
initial_parameters_pc[9, 2] <- dat0422b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[9,4] <- dat0422b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0422b_pc <- dat0422b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0422b_pc <- dat0422b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[9,3] <- (initial_parameters_pc[9,4] - initial_parameters_pc[9,2])/(max_day_0422b_pc - min_day_0422b_pc)

##----0522A: APC----

initial_parameters_apc[10, 1] <- "0522_A"
initial_parameters_apc[10, 5] <- "apc"

#n0 
initial_parameters_apc[10, 2] <- dat0522a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[10,4] <- dat0522a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0522a_apc <- dat0522a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0522a_apc <- dat0522a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[10,3] <- (initial_parameters_apc[10,4] - initial_parameters_apc[10,2])/(max_day_0522a_apc - min_day_0522a_apc)

##----0522A: PC----

initial_parameters_pc[10, 1] <- "0522_A"
initial_parameters_pc[10, 5] <- "pc"

#n0 
initial_parameters_pc[10, 2] <- dat0522a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[10,4] <- dat0522a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0522a_pc <- dat0522a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0522a_pc <- dat0522a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[10,3] <- (initial_parameters_pc[10,4] - initial_parameters_pc[10,2])/(max_day_0522a_pc - min_day_0522a_pc)

##----0522B: APC----

initial_parameters_apc[11, 1] <- "0522_B"
initial_parameters_apc[11, 5] <- "apc"

#n0 
initial_parameters_apc[11, 2] <- dat0522b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[11,4] <- dat0522b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0522b_apc <- dat0522b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0522b_apc <- dat0522b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[11,3] <- (initial_parameters_apc[11,4] - initial_parameters_apc[11,2])/(max_day_0522b_apc - min_day_0522b_apc)

##----0522B: PC----

initial_parameters_pc[11, 1] <- "0522_B"
initial_parameters_pc[11, 5] <- "pc"

#n0 
initial_parameters_pc[11, 2] <- dat0522b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[11,4] <- dat0522b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0522b_pc <- dat0522b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0522b_pc <- dat0522b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[11,3] <- (initial_parameters_pc[11,4] - initial_parameters_pc[11,2])/(max_day_0522b_pc - min_day_0522b_pc)

##----0522C: APC----

initial_parameters_apc[12, 1] <- "0522_C"
initial_parameters_apc[12, 5] <- "apc"

#n0 
initial_parameters_apc[12, 2] <- dat0522c_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[12,4] <- dat0522c_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0522c_apc <- dat0522c_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0522c_apc <- dat0522c_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[12,3] <- (initial_parameters_apc[12,4] - initial_parameters_apc[12,2])/(max_day_0522c_apc - min_day_0522c_apc)

##----0522C: PC----

initial_parameters_pc[12, 1] <- "0522_C"
initial_parameters_pc[12, 5] <- "pc"

#n0 
initial_parameters_pc[12, 2] <- dat0522c_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[12,4] <- dat0522c_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0522c_pc <- dat0522c_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0522c_pc <- dat0522c_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[12,3] <- (initial_parameters_pc[12,4] - initial_parameters_pc[12,2])/(max_day_0522c_pc - min_day_0522c_pc)

##----0622A: APC----

initial_parameters_apc[13, 1] <- "0622_A"
initial_parameters_apc[13, 5] <- "apc"

#n0 
initial_parameters_apc[13, 2] <- dat0622a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[13,4] <- dat0622a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0622a_apc <- dat0622a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0622a_apc <- dat0622a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[13,3] <- (initial_parameters_apc[13,4] - initial_parameters_apc[13,2])/(max_day_0622a_apc - min_day_0622a_apc)

##----0622A: PC----

initial_parameters_pc[13, 1] <- "0622_A"
initial_parameters_pc[13, 5] <- "pc"

#n0 
initial_parameters_pc[13, 2] <- dat0622a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[13,4] <- dat0622a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0622a_pc <- dat0622a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0622a_pc <- dat0622a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[13,3] <- (initial_parameters_pc[13,4] - initial_parameters_pc[13,2])/(max_day_0622a_pc - min_day_0622a_pc)

##----0622B: APC----

initial_parameters_apc[14, 1] <- "0622_B"
initial_parameters_apc[14, 5] <- "apc"

#n0 
initial_parameters_apc[14, 2] <- dat0622b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[14,4] <- dat0622b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0622b_apc <- dat0622b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0622b_apc <- dat0622b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()


#mumax
initial_parameters_apc[14,3] <- (initial_parameters_apc[14,4] - initial_parameters_apc[14,2])/(max_day_0622b_apc - min_day_0622b_apc)

##----0622B: PC----

initial_parameters_pc[14, 1] <- "0622_B"
initial_parameters_pc[14, 5] <- "pc"

#n0 
initial_parameters_pc[14, 2] <- dat0622b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[14,4] <- dat0622b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0622b_pc <- dat0622b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0622b_pc <- dat0622b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[14,3] <- (initial_parameters_pc[14,4] - initial_parameters_pc[14,2])/(max_day_0622b_pc - min_day_0622b_pc)

##----0722A: APC----

initial_parameters_apc[15, 1] <- "0722_A"
initial_parameters_apc[15, 5] <- "apc"

#n0 
initial_parameters_apc[15, 2] <- dat0722a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[15,4] <- dat0722a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0722a_apc <- dat0722a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0722a_apc <- dat0722a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[15,3] <- (initial_parameters_apc[15,4] - initial_parameters_apc[15,2])/(max_day_0722a_apc - min_day_0722a_apc)

##----0722A: PC----

initial_parameters_pc[15, 1] <- "0722_A"
initial_parameters_pc[15, 5] <- "pc"

#n0 
initial_parameters_pc[15, 2] <- dat0722a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[15,4] <- dat0722a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0722a_pc <- dat0722a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0722a_pc <- dat0722a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[15,3] <- (initial_parameters_pc[15,4] - initial_parameters_pc[15,2])/(max_day_0722a_pc - min_day_0722a_pc)

##----0822A: APC----

initial_parameters_apc[16, 1] <- "0822_A"
initial_parameters_apc[16, 5] <- "apc"

#n0 
initial_parameters_apc[16, 2] <- dat0822a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[16,4] <- dat0822a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0822a_apc <- dat0822a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0822a_apc <- dat0822a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[16,3] <- (initial_parameters_apc[16,4] - initial_parameters_apc[16,2])/(max_day_0822a_apc - min_day_0822a_apc)

##----0822A: PC----

initial_parameters_pc[16, 1] <- "0822_A"
initial_parameters_pc[16, 5] <- "pc"

#n0 
initial_parameters_pc[16, 2] <- dat0822a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[16,4] <- dat0822a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0822a_pc <- dat0822a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0822a_pc <- dat0822a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[16,3] <- (initial_parameters_pc[16,4] - initial_parameters_pc[16,2])/(max_day_0822a_pc - min_day_0822a_pc)

##----0822B: APC----

initial_parameters_apc[17, 1] <- "0822_B"
initial_parameters_apc[17, 5] <- "apc"

#n0 
initial_parameters_apc[17, 2] <- dat0822b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[17,4] <- dat0822b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0822b_apc <- dat0822b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0822b_apc <- dat0822b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[17,3] <- (initial_parameters_apc[17,4] - initial_parameters_apc[17,2])/(max_day_0822b_apc - min_day_0822b_apc)

##----0822B: PC----

initial_parameters_pc[17, 1] <- "0822_B"
initial_parameters_pc[17, 5] <- "pc"

#n0 
initial_parameters_pc[17, 2] <- dat0822b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[17,4] <- dat0822b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0822b_pc <- dat0822b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0822b_pc <- dat0822b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[17,3] <- (initial_parameters_pc[17,4] - initial_parameters_pc[17,2])/(max_day_0822b_pc - min_day_0822b_pc)

##----0922A: APC----

initial_parameters_apc[18, 1] <- "0922_A"
initial_parameters_apc[18, 5] <- "apc"

#n0 
initial_parameters_apc[18, 2] <- dat0922a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[18,4] <- dat0922a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0922a_apc <- dat0922a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0922a_apc <- dat0922a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[18,3] <- (initial_parameters_apc[18,4] - initial_parameters_apc[18,2])/(max_day_0922a_apc - min_day_0922a_apc)

##----0922A: PC----

initial_parameters_pc[18, 1] <- "0922_A"
initial_parameters_pc[18, 5] <- "pc"

#n0 
initial_parameters_pc[18, 2] <- dat0922a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[18,4] <- dat0922a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0922a_pc <- dat0922a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0922a_pc <- dat0922a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[18,3] <- (initial_parameters_pc[18,4] - initial_parameters_pc[18,2])/(max_day_0922a_pc - min_day_0922a_pc)

##----0922B: APC----

initial_parameters_apc[19, 1] <- "0922_B"
initial_parameters_apc[19, 5] <- "apc"

#n0 
initial_parameters_apc[19, 2] <- dat0922b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[19,4] <- dat0922b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0922b_apc <- dat0922b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0922b_apc <- dat0922b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[19,3] <- (initial_parameters_apc[19,4] - initial_parameters_apc[19,2])/(max_day_0922b_apc - min_day_0922b_apc)

##----0922B: PC----

initial_parameters_pc[19, 1] <- "0922_B"
initial_parameters_pc[19, 5] <- "pc"

#n0 
initial_parameters_pc[19, 2] <- dat0922b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[19,4] <- dat0922b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_0922b_pc <- dat0922b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_0922b_pc <- dat0922b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[19,3] <- (initial_parameters_pc[19,4] - initial_parameters_pc[19,2])/(max_day_0922b_pc - min_day_0922b_pc)

##----1022A: APC----

initial_parameters_apc[20, 1] <- "1022_A"
initial_parameters_apc[20, 5] <- "apc"

#n0 
initial_parameters_apc[20, 2] <- dat1022a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[20,4] <- dat1022a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1022a_apc <- dat1022a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1022a_apc <- dat1022a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[20,3] <- (initial_parameters_apc[20,4] - initial_parameters_apc[20,2])/(max_day_1022a_apc - min_day_1022a_apc)

##----1022A: PC----

initial_parameters_pc[20, 1] <- "1022_A"
initial_parameters_pc[20, 5] <- "pc"

#n0 
initial_parameters_pc[20, 2] <- dat1022a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[20,4] <- dat1022a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1022a_pc <- dat1022a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1022a_pc <- dat1022a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[20,3] <- (initial_parameters_pc[20,4] - initial_parameters_pc[20,2])/(max_day_1022a_pc - min_day_1022a_pc)

##----1022B: APC----

initial_parameters_apc[21, 1] <- "1022_B"
initial_parameters_apc[21, 5] <- "apc"

#n0 
initial_parameters_apc[21, 2] <- dat1022b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[21,4] <- dat1022b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1022b_apc <- dat1022b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1022b_apc <- dat1022b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[21,3] <- (initial_parameters_apc[21,4] - initial_parameters_apc[21,2])/(max_day_1022b_apc - min_day_1022b_apc)

##----1022B: PC----

initial_parameters_pc[21, 1] <- "1022_B"
initial_parameters_pc[21, 5] <- "pc"

#n0 
initial_parameters_pc[21, 2] <- dat1022b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[21,4] <- dat1022b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1022b_pc <- dat1022b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1022b_pc <- dat1022b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[21,3] <- (initial_parameters_pc[21,4] - initial_parameters_pc[21,2])/(max_day_1022b_pc - min_day_1022b_pc)

##----1122A: APC----

initial_parameters_apc[22, 1] <- "1122_A"
initial_parameters_apc[22, 5] <- "apc"

#n0 
initial_parameters_apc[22, 2] <- dat1122a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[22,4] <- dat1122a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1122a_apc <- dat1122a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1122a_apc <- dat1122a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[22,3] <- (initial_parameters_apc[22,4] - initial_parameters_apc[22,2])/(max_day_1122a_apc - min_day_1122a_apc)

##----1122A: PC----

initial_parameters_pc[22, 1] <- "1122_A"
initial_parameters_pc[22, 5] <- "pc"

#n0 
initial_parameters_pc[22, 2] <- dat1122a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[22,4] <- dat1122a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1122a_pc <- dat1122a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1122a_pc <- dat1122a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[22,3] <- (initial_parameters_pc[22,4] - initial_parameters_pc[22,2])/(max_day_1122a_pc - min_day_1122a_pc)

##----1122B: APC----

initial_parameters_apc[23, 1] <- "1122_B"
initial_parameters_apc[23, 5] <- "apc"

#n0 
initial_parameters_apc[23, 2] <- dat1122b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[23,4] <- dat1122b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1122b_apc <- dat1122b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1122b_apc <- dat1122b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[23,3] <- (initial_parameters_apc[23,4] - initial_parameters_apc[23,2])/(max_day_1122b_apc - min_day_1122b_apc)

##----1122B: PC----

initial_parameters_pc[23, 1] <- "1122_B"
initial_parameters_pc[23, 5] <- "pc"

#n0 
initial_parameters_pc[23, 2] <- dat1122b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[23,4] <- dat1122b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1122b_pc <- dat1122b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1122b_pc <- dat1122b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[23,3] <- (initial_parameters_pc[23,4] - initial_parameters_pc[23,2])/(max_day_1122b_pc - min_day_1122b_pc)

##----1222A: APC----

initial_parameters_apc[24, 1] <- "1222_A"
initial_parameters_apc[24, 5] <- "apc"

#n0 
initial_parameters_apc[24, 2] <- dat1222a_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[24,4] <- dat1222a_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1222a_apc <- dat1222a_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1222a_apc <- dat1222a_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[24,3] <- (initial_parameters_apc[24,4] - initial_parameters_apc[24,2])/(max_day_1222a_apc - min_day_1222a_apc)

##----1222A: PC----

initial_parameters_pc[24, 1] <- "1222_A"
initial_parameters_pc[24, 5] <- "pc"

#n0 
initial_parameters_pc[24, 2] <- dat1222a_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[24,4] <- dat1222a_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1222a_pc <- dat1222a_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1222a_pc <- dat1222a_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[24,3] <- (initial_parameters_pc[24,4] - initial_parameters_pc[24,2])/(max_day_1222a_pc - min_day_1222a_pc)

##----1222B: APC----

initial_parameters_apc[25, 1] <- "1222_B"
initial_parameters_apc[25, 5] <- "apc"

#n0 
initial_parameters_apc[25, 2] <- dat1222b_apc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_apc[25,4] <- dat1222b_apc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1222b_apc <- dat1222b_apc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1222b_apc <- dat1222b_apc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_apc[25,3] <- (initial_parameters_apc[25,4] - initial_parameters_apc[25,2])/(max_day_1222b_apc - min_day_1222b_apc)

##----1222B: PC----

initial_parameters_pc[25, 1] <- "1222_B"
initial_parameters_pc[25, 5] <- "pc"

#n0 
initial_parameters_pc[25, 2] <- dat1222b_pc %>%
  summarise(min(log_conc_geom)) 

#nmax 
initial_parameters_pc[25,4] <- dat1222b_pc %>%
  summarise(max(log_conc_geom)) 

#day with minimum concentration
min_day_1222b_pc <- dat1222b_pc %>%
  arrange(log_conc_geom) %>%
  summarize(day[1]) %>%
  as.numeric()

#day with maximum concentration
max_day_1222b_pc <- dat1222b_pc %>%
  arrange(desc(log_conc_geom)) %>%
  summarize(day[1]) %>%
  as.numeric()

#mumax
initial_parameters_pc[25,3] <- (initial_parameters_pc[25,4] - initial_parameters_pc[25,2])/(max_day_1222b_pc - min_day_1222b_pc)

## -------------------------Counts at the end of shelf life---------------------

counts_end_of_shelf_life <- conc %>%
  filter(test != "GN") %>%
  filter(!(sampling_code %in% c("0422_B", "0522_A", "0922_B"))) %>%
  group_by(loc, sampling_code, test) %>%
  arrange(desc(day)) %>%
  summarize(diff = log_conc_geom[1] - log_conc_geom[2]) %>%
  ungroup() 

counts_end_of_shelf_life$diff_less_0.50 <- ifelse(counts_end_of_shelf_life$diff < 0.50,
                                                  "yes", "no")

## -----------------------------Exporting data----------------------------------

initial_parameters_export_apc <- initial_parameters_apc

initial_parameters_export_pc <- initial_parameters_pc

initial_parameters_export <- bind_rows(initial_parameters_export_apc, initial_parameters_export_pc)

write.csv(initial_parameters_export, "data/wrangled/initial_parameter_estimates.csv", row.names = FALSE)

