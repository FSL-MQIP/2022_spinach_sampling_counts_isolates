## ----------------------Title--------------------------------------------------
#   Plotting sequencing data of isolates

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

# Script description: Plotting the sequencing data of isolates collected in the 
# 2022 sampling. All sequence data collected up until 2023_04_06 was used in the analysis


## ----------------------Packages-----------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(reshape2); library(ggplot2)

## ----------------------Reading in Data----------------------------------------

genus_data <- read.csv("data/wrangled/genus_assignments.csv", header = TRUE)
fmt <- read.csv("data/raw/fmt_Sheet.csv", header = TRUE)
location <- read.csv("data/raw/sample_location_of_origin.csv", header = TRUE)

## ----------------------Wrangling----------------------------------------------

#Assigning the previous ID to the isolate
genus_data$isolate <- paste0("FSL ", toupper(genus_data$qseqid))
fmt$FSL <- toupper(fmt$FSL)

data <- merge(genus_data, fmt, by.x = "isolate", by.y = "FSL", all = TRUE)
data <- data %>%
  select(isolate, genus, Previous.ID) %>%
  filter(!is.na(genus))

colnames(data)[3] <- "previous_id"

data <- separate(data, previous_id, 
                 into = c("sampling_code", "sample", "sample_rep", "test", "plate_rep", "isolate_number", "mixed_culture_no"), 
                 sep = "-")

data$sampling_code_sample <- paste0(data$sampling_code, "_", data$sample)
data$sampling_code_sample_test <- paste0(data$sampling_code, "_", data$sample, "_", data$test)

#Removing data from samples that were transported or tested incorrectly

#H sample, 0722_B
#Shelf life sample, 0922_B
#H sample, 1222A APC and GN only
#Shelf life sample, 0722_B D7 samples -> Not received, wrong isolates frozen down as belonging to this sample
sample_error <- c("0722B_H_APC", "0722B_H_PC", "0722B_H_GN", "0922B_D7_APC", 
                  "0922B_D7_PC", "0922B_D7_GN", "0922B_D22_APC", "0922B_D22_GN",
                  "0922B_D22_PC", "1222A_H_APC", "1222A_H_GN")

data_filt <- data %>%
  filter(!(sampling_code_sample_test %in% sample_error)) 

data_removed <- data %>%
  filter(sampling_code_sample_test %in% sample_error)

#Removing 79 isolates that could not be classified at the genus-level
table((data_filt$genus %in% c("Enterobacteriaceae", "Hafniaceae",
                              "Pectobacteriaeceae", "Micrococcaeae",
                              "Enterobacterales", "Erwiniaceae", 
                              "Caryophanaceae", "Yersiniaceae")))

n_genus <- data_filt %>%
  filter(genus %in% c("Enterobacteriaceae", "Hafniaceae",
                      "Pectobacteriaeceae", "Micrococcaeae",
                      "Enterobacterales", "Erwiniaceae",
                      "Caryophanaceae", "Yersiniaceae"))

data_filt <- data_filt %>%
  filter(!(genus %in% c("Enterobacteriaceae", "Hafniaceae",
                        "Pectobacteriaeceae", "Micrococcaeae",
                        "Enterobacterales", "Erwiniaceae",
                        "Caryophanaceae", "Yersiniaceae")))

genera_prop <- data_filt %>%
  group_by(genus) %>%
  summarize(prop = n()/nrow(data_filt)) %>%
  mutate(percent = prop*100)

genera_high_freq <- genera_prop %>%
  filter(percent > 1)

genera_low_freq <- genera_prop %>%
  filter(percent <= 1)

#Genera with frequency <1% of the total dataset are collapsed into an "rare" category

data_filt$genus_high_freq <- ifelse(data_filt$genus %in% genera_high_freq$genus, 
                                    data_filt$genus,
                                    "Rare")


#Changing D22 to D28 for the following datasets: 0422_A, 0422_B, 0522_A, 
#0522_B, 0522_C, 0622_A

pattern <- c("0422A_D22", "0422B_D22", "0522A_D22", "0522B_D22", "0522C_D22", "0622A_D22")

data_filt$sample <- ifelse(data_filt$sampling_code_sample %in% pattern, "D28", data_filt$sample)

#Adding locations to samples
location$sampling_code <- gsub("_", "", location$sampling_code)
loc_h <- filter(location, sample_type != "f")
loc_f <- filter(location, sample_type != "r")

for(i in 1:nrow(data_filt)){
  if(data_filt$sample[i] == "H"){
    index_h <- which(data_filt$sampling_code[i] == loc_h$sampling_code)
    data_filt$loc[i] <- loc_h$location[index_h]
  } else {
    index_f <- which(data_filt$sampling_code[i] == loc_f$sampling_code)
    data_filt$loc[i] <- loc_f$location[index_f]
  }
}

plot_data <- data_filt %>%
  group_by(loc, sample, genus_high_freq) %>% 
  summarise(genera_count = n()) %>%
  mutate(hf_genera_prop = genera_count/sum(genera_count)) %>%
  mutate(bar_labels = cumsum(hf_genera_prop) - 0.5 * hf_genera_prop)

plot_data$sample <- factor(plot_data$sample, levels = c("H", "D7", "D22", "D28"))

tiff("outputs/genera_by_sample_location_.tiff", units="in", width=7, height=5, res=1200)

ggplot(plot_data, aes(x = sample, y = hf_genera_prop, fill = genus_high_freq)) + 
  geom_bar(stat = "identity", color = "black") + coord_flip() +
  geom_text(aes(label = genera_count), position = position_stack(vjust = 0.5), color = "black", size = 1.75) + 
  labs(y = "Proportion of Isolates", fill = "Genus", x = "Sample") + 
  coord_flip()  + facet_grid(rows = vars(loc)) +
  scale_y_continuous(breaks=c(0,0.5,1)) + 
  theme(strip.text.x = element_text(size = 15), 
        strip.text.y = element_text(size = 15), 
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.4,"cm"),
        legend.key.width = unit(0.4,"cm"), 
        legend.text = element_text(size = 5)) + 
  scale_fill_manual(values = c("#AE76A3", "#1965B0","#7BAFDE", 
                               "#4EB265", "#F6C141", "#666666","#F1932D", "#E8601C")) + 
  theme_bw() + theme(legend.title = element_text(size = 6),legend.key.size = unit(0.5, "cm"),
                     legend.key.width = unit(0.5,"cm"), legend.text = element_text(size = 5), 
                     aspect.ratio = 8/24, legend.position = "top")

dev.off()
