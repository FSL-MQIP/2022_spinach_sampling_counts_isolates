## ------------------------------Title------------------------------------------
# Wrangling Raw Data from the 2022 Spinach Sampling

## ------------------------------Description------------------------------------
# Project: CIDA Spinach 

# Script description: Combining the raw data to a single data frame. Ensuring that 
#                     all variables are in of the appropriate type. Wrangling data as necessary 
#                     (e.g., log transformation). 
#                     

## ------------------------------Packages---------------------------------------
# Loading packages
library(tidyverse); library(dplyr); library(stringr); library(lubridate)

## ------------------------------Raw Data---------------------------------------

sample_mass <- read.csv("data/raw/consolidated_sample_mass.csv")
counts <- read.csv("data/raw/consolidated_counts.csv")
day_initial <- read.csv("data/raw/day_initial_shelf_life.csv")
location <- read.csv("data/raw/sample_location_of_origin.csv")
metadata <- read.csv("data/raw/raw_finished_product_consolidated_metadata_2022sampling.csv", header = TRUE)

## ------------------------------Wrangling--------------------------------------

#Defining a function to calculate the dilution factor for each sample
dilution_factor_func <- function(sample_mass, buffer_mass, plate_dilution){
  first_dilution <- sample_mass/(sample_mass + buffer_mass)
  final_plate_dilution <- (1/first_dilution)*(10^(plate_dilution-1))
  return(final_plate_dilution)
}

#Removing the negative control data from the sampling data
sample_mass_nc <- sample_mass %>%
  filter(grepl("NC", Sample_Name)) %>%
  as.data.frame()

sample_mass_filt_nc <- sample_mass %>%
  filter(!grepl("NC", Sample_Name)) %>%
  as.data.frame()

counts_nc <- counts %>%
  filter(grepl("NC", Sample_Name)) %>%
  as.data.frame()

counts_filt_nc <- counts %>%
  filter(!grepl("NC", Sample_Name)) %>%
  as.data.frame()

#Calculating the sample mass, by subtracting the mass of the empty Whirlpak bag 
#from the total bag mass (Whirlpak bag + spinach)
sample_mass_filt_nc$Sample_Mass_g <- ifelse(is.na(sample_mass_filt_nc$Sample_Mass_g), 
                                            sample_mass_filt_nc$Total_Bag_Mass - sample_mass_filt_nc$Bag_Mass, 
                                            sample_mass_filt_nc$Sample_Mass_g)

#The sample mass was not taken down for 0322B D7. Assuming it is 25g, based on 
#the protocol used to test samples. 
sample_mass_filt_nc$Sample_Mass_g <-ifelse(((sample_mass_filt_nc$Sampling_code == "0322_B") &
                                            (sample_mass_filt_nc$Sample_Name %in% c("D7_1", "D7_2", "D7_3"))), 
                                           25, sample_mass_filt_nc$Sample_Mass_g)

#Calculating the dilution factor for each sample 
counts_filt_nc$dilution_factor <- vector(mode = "logical", length = nrow(counts_filt_nc))


for(i in 1:nrow(counts_filt_nc)){
  sample_index <- which(sample_mass_filt_nc$Sampling_code == counts_filt_nc$Sampling_code[i] & 
                          sample_mass_filt_nc$Sample_Name == counts_filt_nc$Sample_Name[i])
  
  counts_filt_nc$dilution_factor[i] <- dilution_factor_func(sample_mass_filt_nc$Sample_Mass_g[sample_index], 
                                                            sample_mass_filt_nc$Buffer_Volume_mL[sample_index], 
                                                            as.numeric(str_sub(counts_filt_nc$Dilution[i], 5, 5)))
}


#Checking whether each plate replicate was enumerated from the sample dilution, 
# by dividing the dilution factor for each plate replicate by the dilution factor for 
# the other plate replicate. If the dilution factors are the same, then the plate 
#replicates were from the same dilution
#All values were 1, so every pair of enumerated plate replicates were plated with the
#same dilution

plating_dilution_check <- counts_filt_nc %>%
  group_by(Sampling_code, Sample_Name, Test) %>% 
  summarise(df_div = dilution_factor[1]/dilution_factor[2])

table(plating_dilution_check$df_div == 1)

#Calculating the sample concentration, by average the plate counts and multiplying 
#by the dilution factor
counts_final <- counts_filt_nc %>%
  group_by(Sampling_code, Sample_Name, Test) %>%
  mutate(ave_count = mean(Count)) %>%
  ungroup() %>%
  distinct(Sampling_code, Sample_Name, Test, ave_count, dilution_factor) %>%
  mutate(concentration = ave_count * dilution_factor) %>%
  separate(Sample_Name, c("day", "sample_replicate")) 


#GN plates from 0122_B sampling, for H_1, H_2 and H_3, were below the detection limit. 
#Replacing the concentration for these plates with the 25% of the detection limit, and then calculating the log10_concentration 

counts_final$concentration <- ifelse(
  (counts_final$Sampling_code == "0122_B" & counts_final$day == "H" & counts_final$Test == "GN"),
  counts_final$dilution_factor*0.25,
  counts_final$concentration
)

#Calculating the log10-transformed concentration
counts_final <- counts_final %>%
  mutate(log_conc = log10(concentration))

#Replacing day initial with the day in shelf life for all samplings
for(i in 1:nrow(counts_final)){
  index <- which(counts_final$Sampling_code[i] == day_initial$Sampling_Code)
  if(counts_final$day[i] == "DI"){
    counts_final$day[i] <- day_initial$Day[index]
  } else{
    counts_final$day[i] <- counts_final$day[i]
  }
}

#For finished product samples, converting the variable elements (i.e., change D22 to 22)
counts_final$day <- ifelse(counts_final$day != "H", 
                           substr(counts_final$day, 2, nchar(counts_final$day)), counts_final$day)

#Correcting day of testing for lots tested every 7 days, instead of every 5 days
test_7day <- c("0422_A", "0422_B", "0522_A", "0522_B", "0522_C", "0622_A")

counts_final <- counts_final %>%
  mutate(day = case_when((Sampling_code %in% test_7day & day == "12") ~ "14",
                         (Sampling_code %in% test_7day & day == "17") ~ "21",
                         (Sampling_code %in% test_7day & day == "22") ~ "28",
                         TRUE ~ day))

#Removing observations with transit delays or lab errors
#Creating a unique variable for each observation, to facilitate removing observations
counts_final$uniq <- paste0(counts_final$Sampling_code, "_", counts_final$day, "_", counts_final$sample_replicate)

#Removing 0722_B harvest data, since the samples were delayed in 
#transit to Ithaca for >4 days 
samples_to_remove_all <- grep(paste("0722_B_H")
                              , counts_final$uniq, ignore.case = TRUE)

counts_final <- counts_final[-samples_to_remove_all, ]

#Creating a unique variable for each observation, to facilitate removing observations
counts_final$uniq <- paste0(counts_final$Sampling_code, "_", counts_final$day, "_", counts_final$sample_replicate, "_", counts_final$Test)


#Removing 1222_A H APC and GN plates, since they were incubated at 30˚C instead of 35˚C 
samples_to_remove_apcgn <- grep(paste("1222_A_H_1_APC", "1222_A_H_2_APC", "1222_A_H_3_APC", 
                                      "1222_A_H_1_GN", "1222_A_H_2_GN", "1222_A_H_3_GN", 
                                      sep = "|"), counts_final$uniq, ignore.case = TRUE)

counts_final <- counts_final[-samples_to_remove_apcgn, ]

#Adding the location of origin to each sample
raw_prod_orig <- filter(location, sample_type == "r")
finished_prod_orig <- filter(location, sample_type == "f")

counts_final$loc <- vector(mode = "logical", length = nrow(counts_final))

for(i in 1:nrow(counts_final)){
  if(counts_final$day[i] == "H"){
    index <- which(counts_final$Sampling_code[i] == raw_prod_orig$sampling_code)
    counts_final$loc[i] <- raw_prod_orig$location[index]
  } else {
    index <- which(counts_final$Sampling_code[i] == finished_prod_orig$sampling_code)
    counts_final$loc[i] <- finished_prod_orig$location[index]
  }
}

#Replacing the median time of harvest for samples that did not have harvest time recorded
metadata$date_of_harvest_yyyy_mm_dd <- gsub("_", "-", metadata$date_of_harvest_yyyy_mm_dd)

#For 1222B, replacing the date of harvest with 2022_12_09 (harvest occurred over 2022_12_08 and 2022_12_09)

metadata$date_of_harvest_yyyy_mm_dd[metadata$sampling_code == "1222_B"] <- "2022-12-09"

metadata$date_time_harvest_local <- paste0(metadata$date_of_harvest_yyyy_mm_dd, " ", metadata$time_harvest_started_hh.mm)

metadata$date_time_harvest_local <- parse_date_time(metadata$date_time_harvest_local, "YmdHM")

metadata$date_time_harvest_local <- round_date(metadata$date_time_harvest_local, unit = "hour")

harvest_time_not_recorded <- metadata %>%
  filter(time_harvest_started_hh.mm == "not_collected") %>%
  dplyr::select(sampling_code)

independent_samples <- metadata %>%
  filter(paired == "yes") %>%
  group_by(sampling_code) %>%
  slice(1)

independent_samples_2 <- bind_rows(independent_samples, filter(metadata, paired == "no"))

median_time_of_harvest <- independent_samples_2 %>%
  filter(!is.na(date_time_harvest_local)) %>%
  mutate(harvest_hour = as.numeric(gsub(":", "", substr(time_harvest_started_hh.mm, 1, 2))), 
         harvest_minutes = as.numeric(substr(time_harvest_started_hh.mm, 4, 5))) %>%
  arrange(harvest_hour, harvest_minutes) %>%
  ungroup() %>%
  dplyr::select(time_harvest_started_hh.mm) %>%
  slice(median(1:n()))

metadata$time_harvest_started_hh.mm <- ifelse(metadata$sampling_code %in% harvest_time_not_recorded$sampling_code,
                                                       median_time_of_harvest,
                                              metadata$time_harvest_started_hh.mm)

sample_list <- independent_samples_2 %>%
  select(sampling_code, sample_type, paired)

rm(median_time_of_harvest, independent_samples, independent_samples_2, harvest_time_not_recorded)


#FL: Setting timezone to America, NY: Jennings, FL is in the same timezone as NYC, NY
#AZ: Setting timezone to Pheonix, which is in MST (same as Yuma)
#CA: Setting the timezone to Los Angeles, for samples from CA and the 1222_A sample from AZ (from the Imperial Valley)
metadata$date_time_harvest_local <- paste0(metadata$date_of_harvest_yyyy_mm_dd, " ", metadata$time_harvest_started_hh.mm)
metadata$date_time_harvest_az <- ymd_hm(metadata$date_time_harvest_local, tz = "America/Phoenix")
metadata$date_time_harvest_ca <- ymd_hm(metadata$date_time_harvest_local, tz = "America/Los_Angeles")
metadata$date_time_harvest_fl <- ymd_hm(metadata$date_time_harvest_local, tz = "America/New_York")

metadata$date_time_harvest_utc <- ymd_hm(metadata$date_time_harvest_local, tz = "UTC")

for(i in 1:nrow(metadata)){
  loc <- metadata$location[i]
  if(loc == "FL"){
    metadata$date_time_harvest_utc[i] <- with_tz(metadata$date_time_harvest_fl[i], tzone = "UTC")
  } else if(loc == "AZ" & metadata$sampling_code[i] != "1222_A"){
    metadata$date_time_harvest_utc[i] <- with_tz(metadata$date_time_harvest_az[i], tzone = "UTC")
  } else if(loc == "AZ" & metadata$sampling_code[i] == "1222_A"){
    metadata$date_time_harvest_utc[i] <- with_tz(metadata$date_time_harvest_ca[i], tzone = "UTC")
  }else if(loc == "CA"){
    metadata$date_time_harvest_utc[i] <- with_tz(metadata$date_time_harvest_ca[i], tzone = "UTC")
  } else{
    metadata$date_time_harvest_local[i] <-  NA
  }
}

metadata <- metadata %>%
  select(-c("date_time_harvest_az", "date_time_harvest_fl", "date_time_harvest_ca"))

#Changing the timezone of testing to UTC
sample_mass_filt_nc$Date_Time_of_testing_utc <- ymd_hm(sample_mass_filt_nc$Date_Time_of_testing, tz = "America/New_York")

#Separating the harvest data from the shelf life data
harvest_arith <- as.data.frame(filter(counts_final, day == "H"))
shelf_life_arith <- as.data.frame(filter(counts_final, day != "H"))

#Assigning time since harvest to each observation 

harvest_arith$hr_harvest_testing <- vector(mode = "logical", length = nrow(harvest_arith))
harvest_arith$date_of_harvest <- vector(mode = "logical", length = nrow(harvest_arith))

for(i in 1:nrow(harvest_arith)){
  index_metadata <- which(harvest_arith$Sampling_code[i] == metadata$sampling_code & metadata$sample_type == "r")
  index_testing <- which(harvest_arith$Sampling_code[i] == sample_mass_filt_nc$Sampling_code & "H_1" == sample_mass_filt_nc$Sample_Name)
  harvest_arith$hr_harvest_testing[i] <- difftime(sample_mass_filt_nc$Date_Time_of_testing_utc[index_testing], 
                                                  metadata$date_time_harvest_utc[index_metadata], unit = "hours")
  harvest_arith$date_of_harvest[i] <- metadata$date_of_harvest_yyyy_mm_dd[index_metadata]
}

#Calculating days since last irrigation

metadata$date_of_harvest_yyyy_mm_dd <- as.Date(metadata$date_of_harvest_yyyy_mm_dd)
metadata$date_of_last_irrigation_yyyy_mm_dd <- gsub("_", "-", metadata$date_of_last_irrigation_yyyy_mm_dd)
metadata$date_of_last_irrigation_yyyy_mm_dd <- ifelse(metadata$date_of_last_irrigation_yyyy_mm_dd  == "not-collected", 
                                                      NA, metadata$date_of_last_irrigation_yyyy_mm_dd)
metadata$date_of_last_irrigation_yyyy_mm_dd <- as.Date(metadata$date_of_last_irrigation_yyyy_mm_dd)

for(i in 1:nrow(harvest_arith)){
  index_metadata <- which(harvest_arith$Sampling_code[i] == metadata$sampling_code & metadata$sample_type == "r")
  harvest_arith$days_since_last_irrigation[i] <- difftime(as.Date(metadata$date_of_harvest_yyyy_mm_dd[index_metadata]), 
                                                          as.Date(metadata$date_of_last_irrigation_yyyy_mm_dd[index_metadata]), unit = "days")
}

for(i in 1:nrow(shelf_life_arith)){
  index_metadata <- which(shelf_life_arith$Sampling_code[i] == metadata$sampling_code & metadata$sample_type == "f")
  shelf_life_arith$days_since_last_irrigation[i] <- difftime(as.Date(metadata$date_of_harvest_yyyy_mm_dd[index_metadata]), 
                                                          as.Date(metadata$date_of_last_irrigation_yyyy_mm_dd[index_metadata]), unit = "days")
}

#For 0122_A, time of packaging was a duration (18:00 - 18:45). Will replace time of packaging with 18:45
metadata$time_of_packaging[metadata$sampling_code == "0122_A" & metadata$sample_type == "f"] <- "18:45"

#Calculating the time from harvest until packaging and harvest until arrival at the processing plant, for finished product samples

metadata$date_time_packaging_local <- paste0(metadata$date_of_packaging, " ", metadata$time_of_packaging)

metadata$date_time_packaging_local <- parse_date_time(metadata$date_time_packaging_local, "YmdHM", tz = "America/New_York")

metadata$date_time_packaging_local <- round_date(metadata$date_time_packaging_local, "hours")

metadata$date_time_packaging_local <- as_datetime(metadata$date_time_packaging_local)

metadata$date_time_packaging_local <- gsub("UTC", "", metadata$date_time_packaging_local)

metadata$date_time_packaging_utc <- with_tz(metadata$date_time_packaging_local, "UTC")

metadata$date_time_arrival_local <- paste0(gsub("_", "-", metadata$date_of_arrival_yyyy_mm_dd), 
                                           " ", metadata$time_of_arrival_hh.mm)

metadata$date_time_arrival_local <- parse_date_time(metadata$date_time_arrival_local, "YmdHM", tz = "America/New_York")

metadata$date_time_arrival_utc <- with_tz(metadata$date_time_arrival_local, "UTC")

shelf_life_arith$hr_harvest_packaging <- vector(mode = "logical", length = nrow(shelf_life_arith))
shelf_life_arith$hr_harvest_arrival <- vector(mode = "logical", length = nrow(shelf_life_arith))
shelf_life_arith$hr_arrival_packaging <- vector(mode = "logical", length = nrow(shelf_life_arith))

for(i in 1:nrow(shelf_life_arith)){
  sampling_code <- shelf_life_arith$Sampling_code[i]
  index_metadata <- which(sampling_code == metadata$sampling_code & metadata$sample_type == "f")
  shelf_life_arith$hr_harvest_packaging[i] <- difftime(metadata$date_time_packaging_utc[index_metadata], 
                                                       metadata$date_time_harvest_utc[index_metadata], unit = "hours")
  shelf_life_arith$hr_harvest_arrival[i] <- difftime(metadata$date_time_arrival_utc[index_metadata],
                                                     metadata$date_time_harvest_utc[index_metadata], unit = "hours")
  shelf_life_arith$hr_arrival_packaging[i] <- difftime(metadata$date_time_packaging_utc[index_metadata], 
                                                       metadata$date_time_arrival_utc[index_metadata], unit = "hours")
}

#Calculating the geometric mean of each sample
harvest_geom <- harvest_arith %>%
  group_by(Sampling_code, day, Test) %>%
  mutate(log_conc_geom = mean(log_conc)) %>%
  distinct(Sampling_code, day, Test, log_conc_geom, loc, hr_harvest_testing, date_of_harvest, days_since_last_irrigation)

shelf_life_geom <- shelf_life_arith %>%
  group_by(Sampling_code, day, Test) %>%
  mutate(log_conc_geom = mean(log_conc)) %>%
  distinct(Sampling_code, day, Test, log_conc_geom, loc, hr_harvest_packaging, hr_harvest_arrival, hr_arrival_packaging, days_since_last_irrigation)

## ----------------------Preparing data for export------------------------------

sample_mass_nc_wrangled <- sample_mass_nc %>%
  dplyr::select(-c("Date_Time_of_arrival", "Date_Time_of_testing", "Date_Time_of_incubation"))

colnames(sample_mass_nc_wrangled) <- tolower(colnames(sample_mass_nc_wrangled))

counts_nc_wrangled <- counts_nc %>%
  dplyr::select(-c("Date_Time_Enumerated", "Dilution", "Photo_on_labserver_date.", "Comments"))

colnames(counts_nc_wrangled) <- tolower(colnames(counts_nc_wrangled))

harvest_arith <- harvest_arith %>%
  dplyr::select(-c("dilution_factor", "ave_count", "concentration", "uniq"))
colnames(harvest_arith) <- tolower(colnames(harvest_arith))

shelf_life_arith <- shelf_life_arith %>%
  dplyr::select(-c("dilution_factor", "ave_count", "concentration", "uniq"))
colnames(shelf_life_arith) <- tolower(colnames(shelf_life_arith))
  
colnames(harvest_geom) <- tolower(colnames(harvest_geom))
colnames(shelf_life_geom) <- tolower(colnames(shelf_life_geom))

wrangled_metadata <- data.frame(lapply(metadata, as.character), stringsAsFactors = FALSE)

## -----------------------------Exporting data----------------------------------

write.csv(harvest_arith, "data/wrangled/harvest_conc_by_sample_replicate.csv", row.names = FALSE)
write.csv(harvest_geom, "data/wrangled/harvest_conc_by_testing_day.csv", row.names = FALSE)
write.csv(shelf_life_arith, "data/wrangled/shelf_life_conc_by_sample_replicate.csv", row.names = FALSE)
write.csv(shelf_life_geom, "data/wrangled/shelf_life_conc_by_testing_day.csv", row.names = FALSE)
write.csv(sample_mass_nc_wrangled, "data/wrangled/sample_mass_negative_control.csv", row.names = FALSE)
write.csv(counts_nc_wrangled, "data/wrangled/counts_negative_control.csv", row.names = FALSE)
write.csv(wrangled_metadata, "data/wrangled/wrangled_metadata.csv", row.names = FALSE)

