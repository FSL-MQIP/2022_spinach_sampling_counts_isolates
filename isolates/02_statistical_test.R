## ----------------------Title--------------------------------------------------
#   Statistical tests using data gathered from isolates

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

# Script description: Statistical tests using data gathered from isolates


## ----------------------Packages-----------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(reshape2); library(vegan); library(ggpubr); 
library(phyloseq); library(microbiomeMarker); 

## ----------------------Reading in Data----------------------------------------

genus_data <- read.csv("data/wrangled/genus_assignments.csv", header = TRUE)
fmt <- read.csv("data/raw/fmt_sheet.csv", header = TRUE)
location <- read.csv("data/raw/sample_location_of_origin.csv", header = TRUE)
tax_key <- read.csv("data/raw/tax_key.csv", header =TRUE)
  
## ----------------------Wrangling----------------------------------------------
#Assigning the previous ID to the isolate
genus_data$isolate <- paste0("FSL ", toupper(genus_data$qseqid))
fmt$FSL <- toupper(fmt$FSL)

data <- merge(genus_data, fmt, by.x = "isolate", by.y = "FSL", all = TRUE)

#Removing the 171 isolates that could not be sequenced, 
#because they could not be frozen down or failed sequencing

table(is.na(data$genus))

failed_seq_not_frozen <- data %>%
  select(isolate, genus, Previous.ID) %>%
  filter(is.na(genus))

data <- data %>%
  select(isolate, genus, Previous.ID) %>%
  filter(!is.na(genus))

colnames(data)[3] <- "previous_id"

data <- separate(data, previous_id, 
                 into = c("sampling_code", "sample", "sample_rep", "test", "plate_rep", "isolate_number", "mixed_culture_no"), 
                 sep = "-")

data$sampling_code_sample <- paste0(data$sampling_code, "_", data$sample)

#Changing D22 to D28 for the following datasets: 0522_B, 0522_C, 0622_A

pattern <- c("0522B_D22", "0522C_D22", "0622A_D22")

data$sample <- ifelse(data$sampling_code_sample %in% pattern, "D28", data$sample)

data$sampling_code_sample <- paste0(data$sampling_code, "_", data$sample)

data$sampling_code_sample_test <- paste0(data$sampling_code, "_", data$sample, "_", data$test)

#Adding locations to samples
location$sampling_code <- gsub("_", "", location$sampling_code)
loc_h <- filter(location, sample_type != "f")
loc_f <- filter(location, sample_type != "r")

for(i in 1:nrow(data)){
  if(data$sample[i] == "H"){
    index_h <- which(data$sampling_code[i] == loc_h$sampling_code)
    data$loc[i] <- loc_h$location[index_h]
  } else {
    index_f <- which(data$sampling_code[i] == loc_f$sampling_code)
    data$loc[i] <- loc_f$location[index_f]
  }
}


#Removing data from samples that were transported or tested incorrectly

#H sample, 0722_B
#Shelf life sample, 0922_B
#H sample, 1222A APC and GN only
#Shelf life sample, 0722_B D7 samples -> Not received, wrong isolates frozen down as belonging to this sample
sample_error <- c("0722B_H_APC", "0722B_H_PC", "0722B_H_GN", "0922B_D7_APC", 
                  "0922B_D7_PC", "0922B_D7_GN", "0922B_D22_APC", "0922B_D22_GN",
                  "0922B_D22_PC", "1222A_H_APC", "1222A_H_GN",
                  "0722B_D7_APC", "0722B_D7_GN", "0722B_D7_PC")

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
  summarize(prop = n()/nrow(data_filt), 
            count = n()) %>%
  mutate(percent = prop*100)

genera_high_freq <- genera_prop %>%
  filter(percent > 1)

#Genera with frequency <1% of the total dataset are collapsed into an "rare" category

data_filt$genus_high_freq <- ifelse(data_filt$genus %in% genera_high_freq$genus, 
                                    data_filt$genus,
                                    "rare")

data_filt_summary_loc <- data_filt %>% 
  group_by(loc) %>%
  summarize(n = n())

data_filt_summary_samp <- data_filt %>% 
  group_by(sample) %>%
  summarize(n = n())

data_filt_summary_genus <- data_filt %>% 
  group_by(genus_high_freq) %>%
  summarize(n = n())

data_filt_summary_loc_genus <- data_filt %>% 
  group_by(loc, genus_high_freq) %>%
  summarize(n = n())

data_filt_summary_loc_samp_num_isolates <- data_filt %>% 
  group_by(loc, sample, genus_high_freq) %>%
  summarize(n = n())

data_filt_summary_samp_genera <- data_filt %>% 
  group_by(sample, genus_high_freq) %>%
  summarize(n = n())


#Chi-square test: Location

chisq.test(data_filt$loc, data_filt$genus_high_freq, simulate.p.value = TRUE)$expected

#Expected value is greater than 5 in at less 80% of the cells. Will use a fisher's exact test subsequently

fisher.test(data_filt$loc, data_filt$genus_high_freq, simulate.p.value = TRUE)

#Chi-square test: Sample

fisher.test(data_filt$sample, data_filt$genus_high_freq, simulate.p.value = TRUE)


## ----------------------Wrangling----------------------------------------------

#Samples tested every 7 days 
seven_days <- c("0422B", "0522A", "0522B", "0522C", "0622A")
  
data_filt_2 <- data_filt %>%
  filter(!(sampling_code %in% seven_days & sample %in% c("D28"))) %>%
  filter(!(loc %in% c("FL", "GA")))

dist_data <- select(data_filt_2, c("genus_high_freq", "sampling_code_sample", "test"))

genus_by_sample <- melt(dist_data, id.vars = c("sampling_code_sample", "genus_high_freq"))

com_mat <- dcast(genus_by_sample, sampling_code_sample ~ genus_high_freq, fun.aggregate = length)

com_mat_met <- separate(com_mat, sampling_code_sample, into = c("sampling_code", "sample"), sep = "_") %>%
  select(sampling_code, sample)
row.names(com_mat_met) <- 1:nrow(com_mat_met)

for(i in 1:nrow(com_mat_met)){
  if(com_mat_met$sample[i] == "H"){
    index_h <- which(com_mat_met$sampling_code[i] == loc_h$sampling_code)
    com_mat_met$loc[i] <- loc_h$location[index_h]
  } else {
    index_f <- which(com_mat_met$sampling_code[i] == loc_f$sampling_code)
    com_mat_met$loc[i] <- loc_f$location[index_f]
  }
}

com_mat_met$sample <- as.factor(com_mat_met$sample)

row.names(com_mat) <- com_mat[, 1]

bc_mat <- vegdist(select(com_mat, !"sampling_code_sample"), binary = FALSE, method = "bray", diag = TRUE)

pcoa <- ape::pcoa(bc_mat)
df.pcoa <- data.frame(pcoa$vectors)
pcoa$values

plot_data <- cbind(com_mat_met, df.pcoa)

x_lab <- round(pcoa$values["Relative_eig"][[1]][1]*100, 2)
y_lab <- round(pcoa$values["Relative_eig"][[1]][2]*100, 2)
z_lab <- round(pcoa$values["Relative_eig"][[1]][3]*100, 2)

#By location only 
cols <- c("#0077bb", "#ddaa33", "#009988", "#cc3311")

loc_pcoa <- ggplot(plot_data,aes(x=Axis.1,y=Axis.2, color = loc))+ 
  geom_point(size=2) + theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "top", legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(0,0,0,0),  legend.title=element_text(size=8), 
        legend.text=element_text(size=7)) + 
  labs(x=paste0("PCoA1 (",x_lab,"%)"), y=paste0("PCoA2 (",y_lab,"%)")) + 
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = cols) + labs(color = "Growing region") +
  stat_ellipse(type='t', level = 0.95)

#By sample only 
sample_pcoa <- ggplot(plot_data,aes(x=Axis.1,y=Axis.2, color = sample))+ 
  geom_point(size=2) + theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "top", legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),  legend.title=element_text(size=8), 
        legend.text=element_text(size=7)) + 
  labs(x=paste0("PCoA1 (",x_lab,"%)"), y=paste0("PCoA2 (",y_lab,"%)")) + 
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_brewer(palette = "Dark2") + labs(color = "Sample") + 
  stat_ellipse(type='t', level = 0.95)

pcoa_plot <- ggarrange(loc_pcoa + rremove("ylab") + rremove("xlab"), 
                       sample_pcoa + rremove("ylab") + rremove("xlab"), labels = c("A", "B"))


## ----------------------Analysis-----------------------------------------------
set.seed(1)
perm_an <- adonis2(bc_mat ~ sample + loc + loc : sample, 
                   data = com_mat_met, permutations = 999, by = "terms")
perm_an

#Posthoc: LDA, sample and timepoint 

com_mat_row <- dcast(genus_by_sample, genus_high_freq ~ sampling_code_sample,
                     fun.aggregate = length)

com_mat_row$genus_high_freq <- as.factor(com_mat_row$genus_high_freq)
rownames(com_mat_row) <- com_mat_row$genus_high_freq 
com_mat_row <- com_mat_row[, -1]
com_mat_row <- as.matrix(com_mat_row)
class(com_mat_row) <- "numeric"

tax_ply <- tax_key 

tax_ply[6, ] <- "rare"


tax_ply$genus <- as.factor(tax_ply$analysis_genera_lda)
tax_ply$species <- "sp"
tax_ply <- tax_ply %>%
  dplyr::select("domain", "phylum", "class", "order", "family", "genus", "species")
tax_ply <- as.data.frame(tax_ply)
rownames(tax_ply) <- tax_ply[, 6]
colnames(tax_ply) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_ply <- as.matrix(tax_ply)

com_mat_met_phy <- com_mat_met
rownames(com_mat_met_phy) <- paste0(com_mat_met_phy$sampling_code, "_", com_mat_met_phy$sample)
com_mat_met_phy$loc <- as.factor(com_mat_met_phy$loc)
com_mat_met_phy$loc <- relevel(com_mat_met_phy$loc, ref = "CA")

OTU = otu_table(com_mat_row, taxa_are_rows = TRUE)
TAX = tax_table(tax_ply)
SAMP = sample_data(com_mat_met_phy)

phy_obj = phyloseq(OTU, TAX, SAMP)

lefse <- run_lefse(phy_obj, 
                   taxa_rank = "Genus",
                   group = "loc",
                   lda_cutoff = 4, 
                   multigrp_strat = TRUE)
plot_ef_bar(lefse, label_level = 2, max_label_len = 60)

lefse <- run_lefse(phy_obj, 
                   taxa_rank = "Genus",
                   group = "sample", 
                   lda_cutoff = 4,
                   multigrp_strat = TRUE)
plot_ef_bar(lefse, label_level = 2, max_label_len = 60)

