## ------------------------------Title------------------------------------------
# Plotting the sampling data 

## ------------------------------Description------------------------------------
# Project: CIDA Spinach 

# Script description: Plotting the sampling data 

## ------------------------------Packages---------------------------------------
# Loading packages
library(tidyverse); library(dplyr); library(ggplot2); library(RColorBrewer); 
library(ggpubr); library(grid); library(reshape2); library(gridExtra)

## ------------------------------Raw Data---------------------------------------

harvest_arith <- read.csv("data/wrangled/harvest_conc_by_sample_replicate.csv", header = TRUE)
harvest_geom <- read.csv("data/wrangled/harvest_conc_by_testing_day.csv", header = TRUE)

shelf_life_arith <- read.csv("data/wrangled/shelf_life_conc_by_sample_replicate.csv", header = TRUE)
shelf_life_geom <- read.csv("data/wrangled/shelf_life_conc_by_testing_day.csv", header = TRUE)

metadata <- read.csv("data/wrangled/wrangled_metadata.csv", header = TRUE)

param <- read.csv("output/primary_growth_model_parameters.csv", header = TRUE)

## ------------------------------Defining models--------------------------------

#Baranyi without lag
baranyi_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10nmax - log10(1 + (10^(log10nmax - log10n0) - 1) * exp(-mumax * day))
}

#Buchanan without lag
buchanan_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10n0 + (day <= ((log10nmax - log10n0) * log(10) / mumax)) * mumax * day / log(10) + (day > ((log10nmax - log10n0) * log(10) / mumax)) * (log10nmax - log10n0)
}

growth_func <- function(day, log10n0, log10nmax, mumax, model){
  if(model == "buchanan_no_lag"){
    return(buchanan_no_lag(day, log10n0, log10nmax, mumax))
  } else if(model == "baranyi_no_lag"){
    return(baranyi_no_lag(day, log10n0, log10nmax, mumax))
  }
}

## ------------------------------Wrangling--------------------------------------

#Converting days into a numeric vector, for the shelf life dataframes
shelf_life_arith$day <- as.numeric(shelf_life_arith$day)

shelf_life_geom$day <- as.numeric(shelf_life_geom$day)


#Sorting samples by lot, for plotting
shelf_life_arith$sampling_code <- as.factor(shelf_life_arith$sampling_code)

shelf_life_arith$sample_replicate <- as.factor(shelf_life_arith$sample_replicate)

arith_1221a <- filter(shelf_life_arith, sampling_code == "1221_A")

arith_0122a <- filter(shelf_life_arith, sampling_code == "0122_A")
arith_0122b <- filter(shelf_life_arith, sampling_code == "0122_B")

arith_0222a <- filter(shelf_life_arith, sampling_code == "0222_A")
arith_0222b <- filter(shelf_life_arith, sampling_code == "0222_B")

arith_0322a <- filter(shelf_life_arith, sampling_code == "0322_A")
arith_0322b <- filter(shelf_life_arith, sampling_code == "0322_B")

arith_0422a <- filter(shelf_life_arith, sampling_code == "0422_A")
arith_0422b <- filter(shelf_life_arith, sampling_code == "0422_B")

arith_0522a <- filter(shelf_life_arith, sampling_code == "0522_A")
arith_0522b <- filter(shelf_life_arith, sampling_code == "0522_B")
arith_0522c <- filter(shelf_life_arith, sampling_code == "0522_C")


arith_0622a <- filter(shelf_life_arith, sampling_code == "0622_A")
arith_0622b <- filter(shelf_life_arith, sampling_code == "0622_B")

arith_0722a <- filter(shelf_life_arith, sampling_code == "0722_A")

arith_0822a <- filter(shelf_life_arith, sampling_code == "0822_A")
arith_0822b <- filter(shelf_life_arith, sampling_code == "0822_B")

arith_0922a <- filter(shelf_life_arith, sampling_code == "0922_A")
arith_0922b <- filter(shelf_life_arith, sampling_code == "0922_B")

arith_1022a <- filter(shelf_life_arith, sampling_code == "1022_A")
arith_1022b <- filter(shelf_life_arith, sampling_code == "1022_B")

arith_1122a <- filter(shelf_life_arith, sampling_code == "1122_A")
arith_1122b <- filter(shelf_life_arith, sampling_code == "1122_B")

arith_1222a <- filter(shelf_life_arith, sampling_code == "1222_A")
arith_1222b <- filter(shelf_life_arith, sampling_code == "1222_B")

## -----------------------------Plotting----------------------------------------
#Plotting harvest samples by location, for AZ and C

harvest_plot_data <- harvest_geom %>%
  filter(!(sampling_code %in% c("0722_B"))) %>%
  filter(!(sampling_code == "1222_A" & test %in% c("APC", "GN")))

harvest_plot_data$year_month_harvest <- substr(harvest_plot_data$date_of_harvest,1, 7)

cols <- c("#0077bb", "#ddaa33", "#009988", "#cc3311")

#Harvest counts by month of collection
harvest_by_sampling_date <- ggplot(data = harvest_plot_data, mapping = aes(x = year_month_harvest, y = log_conc_geom), color = loc) + 
 geom_jitter(aes(color = loc), width = 0.1, height = 0, size = 0.75) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_color_manual(values = cols) + 
  ylab(expression("Log"[10]*"CFU/g")) + xlab("Month of sampling, YYYY-MM") + scale_fill_manual(values = cols) + facet_wrap(facets = vars(test), nrow = 3) + labs(colour="Growing region")


harvest_by_sampling_date_boxplot <- ggplot(data = filter(harvest_plot_data, loc %in% c("AZ", "CA")), mapping = aes(x = loc, y = log_conc_geom), fill = loc) + 
  geom_boxplot(aes(fill = loc)) + ylab(expression("Log"[10]*"CFU/g")) + xlab("Growing region") + scale_fill_manual(values = cols)  + facet_wrap(facets = vars(test), nrow = 3) + theme_bw() + labs(colour="Growing region")

#tiff("output/harvest_plot.tiff", width = 6, height = 4, unit = "in", res = 1200)
ggarrange(harvest_by_sampling_date, 
          harvest_by_sampling_date_boxplot + rremove("ylab") + rremove("xlab"), labels = c("A", "B"),
          widths = c(5,2), heights = c(1, 3), common.legend = TRUE, legend = "right",
          align = "hv")
#dev.off()

#D7 counts by month of collection

shelf_life_plot_data <- shelf_life_geom %>%
  filter(!(sampling_code %in% c("0922_B")))

shelf_life_plot_data$year_month_harvest <- paste0("20",substr(shelf_life_plot_data$sampling_code, 3, 4), 
                                             "-", substr(shelf_life_plot_data$sampling_code, 1, 2))
                                             
d7_by_sampling_date <- ggplot(data = filter(shelf_life_plot_data, day == "7"), mapping = aes(x = year_month_harvest, y = log_conc_geom), color = loc) + 
  geom_jitter(aes(color = loc), width = 0.1, height = 0, size = 0.75) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_color_manual(values = cols) + 
  ylab(expression("Log"[10]*"CFU/g")) + xlab("Month of sampling, YYYY-MM") + facet_wrap(facets = vars(test), nrow = 3) + labs(colour="Growing region")

d7_by_sampling_date_boxplot <- ggplot(data = filter(shelf_life_plot_data, day == "7", loc %in% c("AZ", "CA")), mapping = aes(x = loc, y = log_conc_geom), fill = loc) + 
  geom_boxplot(aes(fill = loc)) + ylab(expression("Log"[10]*"CFU/g"))  + xlab("Growing region")  + scale_fill_manual(values = cols)  + facet_wrap(facets = vars(test), nrow = 3) + theme_bw() + labs(colour="Growing region")

d7_plots <- ggarrange(d7_by_sampling_date, 
          d7_by_sampling_date_boxplot + rremove("ylab") + rremove("xlab"), labels = c("A", "B"),
          nrow = 1, ncol = 2, align = "hv", widths = c(5,2), heights = c(1, 3),
          common.legend = TRUE, legend = "right")

#tiff("output/d7_plot.tiff", width = 6, height = 4, unit = "in", res = 1200)
d7_plots
#dev.off()

#Plotting shelf life data together 
shelf_life_plot <- ggplot(data = filter(shelf_life_plot_data), mapping = aes(x = day, y = log_conc_geom)) + 
  geom_jitter(aes(color = loc),shape = 16, width = 0.2, size = 0.75) + facet_grid(rows = vars(test)) + 
  scale_color_manual(values = cols) + labs(colour="Growing region") + 
  ylab(expression("Log"[10]*"CFU/g")) + xlab("Day") + theme_bw()

#tiff("output/shelf_life_plot.tiff", width = 6, height = 4, unit = "in", res = 1200)
shelf_life_plot
#dev.off()

## -----------------------------Plotting----------------------------------------

samp <- rep(unique(shelf_life_arith$sampling_code), each = 29)
day <- c(rep(0:28, times = 25))
pred <- vector(mode = "logical", length = length(day))

apc_pred <- bind_cols(samp, day, pred)
colnames(apc_pred) <- c("sampling_code", "day", "predicted_concentration")

pc_pred <- apc_pred

for(i in 1:nrow(apc_pred)){
  if(!(apc_pred$sampling_code[i] %in% c("0422_B", "0522_A"))){
  index <- which(apc_pred$sampling_code[i] == param$sampling_code &
                   "apc" == param$test)
  apc_pred$predicted_concentration[i] <- growth_func(apc_pred$day[i], 
                                                  param$n0[index], 
                                                  param$nmax[index], 
                                                  param$mumax[index], 
                                                  param$model[index])
  } else {
    apc_pred$predicted_concentration[i] <- NA
  }}

for(i in 1:nrow(pc_pred)){
  if(!(pc_pred$sampling_code[i] %in% c("0422_B", "0522_A"))){
    index <- which(pc_pred$sampling_code[i] == param$sampling_code &
                     "pc" == param$test)
    pc_pred$predicted_concentration[i] <- growth_func(pc_pred$day[i], 
                                                       param$n0[index], 
                                                       param$nmax[index], 
                                                       param$mumax[index], 
                                                       param$model[index])
  } else {
    pc_pred$predicted_concentration[i] <- NA
  }}

#APC 

apc_0122a <- ggplot(data = filter(arith_0122a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0122_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0122_A, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0122b <- ggplot(data = filter(arith_0122b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0122_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0122_B, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)


apc_0222a <- ggplot(data = filter(arith_0222a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0222_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0222_A, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0222b <- ggplot(data = filter(arith_0222b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0222_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0222_B, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0322a <- ggplot(data = filter(arith_0322a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0322_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0322_A, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0322b <- ggplot(data = filter(arith_0322b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0322_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0322_B, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0422a <- ggplot(data = filter(arith_0422a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0422_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0422_A, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0522b <- ggplot(data = filter(arith_0522b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0522_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0522_B, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0522c <- ggplot(data = filter(arith_0522c, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0522_C"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0522_C, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0622a <- ggplot(data = filter(arith_0622a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0622_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0622_A, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0622b <- ggplot(data = filter(arith_0622b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0622_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0622_B, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0722a <- ggplot(data = filter(arith_0722a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0722_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0722_A, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0822a <- ggplot(data = filter(arith_0822a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0822_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0822_A, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0822b <- ggplot(data = filter(arith_0822b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0822_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0822_B, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0922a <- ggplot(data = filter(arith_0922a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0922_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0922_A, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_0922b <- ggplot(data = filter(arith_0922b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "0922_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0922_B, APC  (NA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_1022a <- ggplot(data = filter(arith_1022a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "1022_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1022_A, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_1022b <- ggplot(data = filter(arith_1022b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "1022_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1022_B, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_1122a <- ggplot(data = filter(arith_1122a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "1122_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1122_A, APC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_1122b <- ggplot(data = filter(arith_1122b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "1122_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1122_B, APC (GA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_1221a <- ggplot(data = filter(arith_1221a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "1221_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1221_A, APC (FL)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_1222a <- ggplot(data = filter(arith_1222a, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "1222_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1222_A, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_1222b <- ggplot(data = filter(arith_1222b, test == "APC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(apc_pred, sampling_code == "1222_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1222_B, APC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5), 
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10.0)

apc_plot <- ggarrange(apc_1221a + rremove("ylab") + rremove("xlab"), 
                         apc_0122a +rremove("ylab") + rremove("xlab"),
                         apc_0122b +rremove("ylab") + rremove("xlab"), 
                         apc_0222a +rremove("ylab") + rremove("xlab"),
                         apc_0222b +rremove("ylab") + rremove("xlab"), 
                         apc_0322a +rremove("ylab") + rremove("xlab"),
                         apc_0322b +rremove("ylab") + rremove("xlab"), 
                         apc_0422a +rremove("ylab") + rremove("xlab"),
                         apc_0522b +rremove("ylab") + rremove("xlab"), 
                         apc_0522c +rremove("ylab") + rremove("xlab"),
                         apc_0622a +rremove("ylab") + rremove("xlab"), 
                         apc_0622b +rremove("ylab") + rremove("xlab"),
                         apc_0722a +rremove("ylab") + rremove("xlab"), 
                         apc_0822a +rremove("ylab") + rremove("xlab"),
                         apc_0822b +rremove("ylab") + rremove("xlab"), 
                         apc_0922a +rremove("ylab") + rremove("xlab"),
                         apc_0922b +rremove("ylab") + rremove("xlab"), 
                         apc_1022a +rremove("ylab") + rremove("xlab"),
                         apc_1022b +rremove("ylab") + rremove("xlab"), 
                         apc_1122a +rremove("ylab") + rremove("xlab"),
                         apc_1122b +rremove("ylab") + rremove("xlab"), 
                         apc_1222a +rremove("ylab") + rremove("xlab"),
                         apc_1222b +rremove("ylab") + rremove("xlab"),
                         common.legend = TRUE, 
                         ncol = 3, nrow = 9)


pdf("output/apc_primary_model_plot.pdf")

annotate_figure(apc_plot, left = textGrob(expression("Concentration, log"[10]*"CFU/g"), rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                bottom = textGrob("Day", y = 3.5, x = 0.5, gp = gpar(cex = 0.8)))

dev.off()

#PC 

pc_0122a <- ggplot(data = filter(arith_0122a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0122_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0122_A, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0122b <- ggplot(data = filter(arith_0122b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0122_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0122_B, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)


pc_0222a <- ggplot(data = filter(arith_0222a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0222_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0222_A, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0222b <- ggplot(data = filter(arith_0222b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0222_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0222_B, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0322a <- ggplot(data = filter(arith_0322a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0322_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0322_A, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0322b <- ggplot(data = filter(arith_0322b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0322_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0322_B, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0422a <- ggplot(data = filter(arith_0422a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0422_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0422_A, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0522b <- ggplot(data = filter(arith_0522b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0522_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0522_B, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0522c <- ggplot(data = filter(arith_0522c, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0522_C"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0522_C, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0622a <- ggplot(data = filter(arith_0622a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0622_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0622_A, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0622b <- ggplot(data = filter(arith_0622b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0622_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0622_B, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0722a <- ggplot(data = filter(arith_0722a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0722_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0722_A, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0822a <- ggplot(data = filter(arith_0822a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0822_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0822_A, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0822b <- ggplot(data = filter(arith_0822b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0822_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0822_B, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0922a <- ggplot(data = filter(arith_0922a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0922_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0922_A, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_0922b <- ggplot(data = filter(arith_0922b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "0922_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "0922_B, PC (NA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_1022a <- ggplot(data = filter(arith_1022a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "1022_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1022_A, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_1022b <- ggplot(data = filter(arith_1022b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "1022_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1022_B, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_1122a <- ggplot(data = filter(arith_1122a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "1122_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1122_A, PC (CA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_1122b <- ggplot(data = filter(arith_1122b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "1122_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1122_B, PC (GA)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_1221a <- ggplot(data = filter(arith_1221a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "1221_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1221_A, PC (FL)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_1222a <- ggplot(data = filter(arith_1222a, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "1222_A"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1222_A, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_1222b <- ggplot(data = filter(arith_1222b, test == "PC"), 
                    aes(x = day, y = log_conc)) + 
  geom_point(aes(shape = sample_replicate), size = 0.5) + 
  geom_line(data = filter(pc_pred, sampling_code == "1222_B"), 
            aes(x = day, y = predicted_concentration), linewidth = 0.2) + 
  xlab("Day") + ylab(expression("Log"[10]*"CFU/g")) + 
  labs(title = "1222_B, PC (AZ)") + 
  scale_shape_manual(values = c(0, 1, 2)) + 
  labs(shape = "Technical replicate") + 
  theme_bw() +
  theme(plot.title = element_text(vjust = -3.5, size = 5),
        axis.text=element_text(size=5), 
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 4)) +
  ylim(3.8, 10)

pc_plot <- ggarrange(pc_1221a + rremove("ylab") + rremove("xlab"), 
                      pc_0122a +rremove("ylab") + rremove("xlab"),
                      pc_0122b +rremove("ylab") + rremove("xlab"), 
                      pc_0222a +rremove("ylab") + rremove("xlab"),
                      pc_0222b +rremove("ylab") + rremove("xlab"), 
                      pc_0322a +rremove("ylab") + rremove("xlab"),
                      pc_0322b +rremove("ylab") + rremove("xlab"), 
                      pc_0422a +rremove("ylab") + rremove("xlab"),
                      pc_0522b +rremove("ylab") + rremove("xlab"), 
                      pc_0522c +rremove("ylab") + rremove("xlab"),
                      pc_0622a +rremove("ylab") + rremove("xlab"), 
                      pc_0622b +rremove("ylab") + rremove("xlab"),
                      pc_0722a +rremove("ylab") + rremove("xlab"), 
                      pc_0822a +rremove("ylab") + rremove("xlab"),
                      pc_0822b +rremove("ylab") + rremove("xlab"), 
                      pc_0922a +rremove("ylab") + rremove("xlab"),
                      pc_0922b +rremove("ylab") + rremove("xlab"), 
                      pc_1022a +rremove("ylab") + rremove("xlab"),
                      pc_1022b +rremove("ylab") + rremove("xlab"), 
                      pc_1122a +rremove("ylab") + rremove("xlab"),
                      pc_1122b +rremove("ylab") + rremove("xlab"), 
                      pc_1222a +rremove("ylab") + rremove("xlab"),
                      pc_1222b +rremove("ylab") + rremove("xlab"),
                      common.legend = TRUE, 
                      ncol = 3, nrow = 9)


pdf("output/pc_primary_model_plot.pdf")

annotate_figure(pc_plot, left = textGrob(expression("Concentration, log"[10]*"CFU/g"), rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                bottom = textGrob("Day", y = 3.5, x = 0.5, gp = gpar(cex = 0.8)))

dev.off()

