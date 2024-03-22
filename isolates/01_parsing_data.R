## ----------------------Title--------------------------------------------------
#   Parsing BLAST data

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

# Script description: Parsing the BLAST data of the isolates from the 
# 2022 sampling.


## ----------------------Packages-----------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(reshape2)

## ----------------------Reading in Data----------------------------------------

genus <- read.csv("data/raw/blast_results.csv", header = FALSE)

## ----------------------Wrangling----------------------------------------------

colnames(genus) <- c("qseqid", "sseqid", "qlen", "slen", "qstart", "qend", 
                     "length", "mismatch", "gapopen", "gaps", "evalue", "bitscore", 
                     "pident", "nident", "sseq")

#Only retaining BLAST hits where the alignment in the BLAST match was
#at least 99% of the total sequence length 

genus_filt <- genus %>% 
  mutate(align_length = length/qlen) %>%
  filter(align_length > 0.99)
  

#Identify which isolates had a high hit for two or more genera

mult_genera <- genus_filt %>%
  separate(sseqid, into = c("genus", NA, NA, NA)) %>%
  group_by(qseqid, genus) %>%
  arrange(desc(pident)) %>%
  slice(1) %>%
  select(qseqid, genus, pident, nident) %>%
  group_by(qseqid) %>%
  summarize(obs = n()) %>%
  filter(obs > 1)

mult_genera_exp <- genus_filt %>%
  separate(sseqid, into = c("genus", NA, NA, NA)) %>%
  group_by(qseqid, genus) %>%
  arrange(desc(pident)) %>%
  slice(1) %>%
  ungroup() %>%
  select(qseqid, genus, pident, nident) %>%
  filter(qseqid %in% mult_genera$qseqid) %>%
  group_by(qseqid) %>%
  arrange(desc(pident)) %>%
  mutate(diff_pident = pident[1] - pident)%>%
  mutate(diff_nident = nident[1] - nident)%>%
  ungroup() 

isolates_to_check <- mult_genera_exp %>%
  filter(diff_pident <= 0.5) %>%
  filter(duplicated(qseqid))

blast_results_mult <- genus_filt %>%
  filter(qseqid %in% isolates_to_check$qseqid) %>%
  group_by(qseqid) %>%
  arrange(desc(pident)) %>%
  mutate(diff_pident = pident[1] - pident, 
         diff_nident = nident[1] - nident)%>% 
  filter(diff_pident <= 0.5) %>%
  ungroup()

#Classifying isolates at a higher taxonomic level, for isolates that 
#had a hit in a different genera than the top hit that was within 0.55 percentage points
#of the percent identity of the top hit 

family_enterobacteriaceae <- c("s12-0468", "s12-0679", "s12-0762", "s12-0949",
                               "s12-0961", "s12-1477", "s12-1478", "s12-1479",
                               "s12-1480", "s12-1513", "s12-1528", "s12-1529",
                               "s12-1890", "s12-1911", "s12-2001", "s12-2051",
                               "s12-2064", "s12-2066", "s12-2099", "s12-2191",
                               "s12-2193", "s12-2194")

family_hafniaceae <- c("s12-0940", "s12-0983", "s12-1103")

family_pectobacteriaceae <- c("s12-1306")

family_erwiniaceae <- c("s12-0457", "s12-0459", "s12-0469",
                        "s12-0479", "s12-0481", "s12-0673", "s12-1129",
                        "s12-2622")

family_micrococcaceae <- c("s12-2225")

family_caryophanaceae <- c("s12-0958")

family_yersiniaceae <- c("s12-2451")

order_enterobacterales <- c("s12-0339", "s12-0427", "s12-0811", "s12-1173",
                            "s12-1177", "s12-1178", "s12-1186", "s12-1189",
                            "s12-1194", "s12-1195", "s12-1196", "s12-1197",
                            "s12-1198", "s12-1215", "s12-1223", "s12-1317",
                            "s12-1318", "s12-1319", "s12-1320", "s12-1325",
                            "s12-1326", "s12-1327", "s12-1328", "s12-1334",
                            "s12-1335", "s12-1336", "s12-1392", "s12-1395",
                            "s12-1406", "s12-1525", "s12-1526", "s12-1579",
                            "s12-1631", "s12-1690", "s12-1969", "s12-1973",
                            "s12-2030", "s12-2033", "s12-2054", "s12-2405",
                            "s12-2411", "s12-2464", "s12-2498", "s12-2551",
                            "s12-2652", "s12-2663", "s12-2670")

genus_final <- genus_filt %>%
  group_by(qseqid) %>%
  arrange(desc(pident)) %>%
  slice(1) %>%
  select(qseqid, sseqid) %>%
  separate(sseqid, into = c("genus", NA, NA, NA))

genus_final <- genus_final %>%
  mutate(genus = case_when(
    qseqid %in% family_enterobacteriaceae ~ "Enterobacteriaceae",
    qseqid %in% family_hafniaceae ~ "Hafniaceae",
    qseqid %in% family_pectobacteriaceae ~ "Pectobacteriaeceae",
    qseqid %in% family_micrococcaceae ~ "Micrococcaeae",
    qseqid %in% order_enterobacterales ~ "Enterobacterales",
    qseqid %in% family_erwiniaceae ~ "Erwiniaceae",
    qseqid %in% family_caryophanaceae ~ "Caryophanaceae",
    qseqid %in% family_yersiniaceae ~ "Yersiniaceae",
    .default = genus
  ))



## ----------------------Export-------------------------------------------------


write.csv(genus_final, "data/wrangled/genus_assignments.csv", row.names = FALSE)

