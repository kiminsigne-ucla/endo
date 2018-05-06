#This script takes Barcode RNA and DNA counts, 
# merges them to their respective promoters, 
# and calculates the strength of each promoter

#Install required packages
# install.packages("dplyr")
# install.packages("ggplot2")

library(dplyr)
library(tidyr)
options(stringsAsFactors = F)

#SET TO WORKING DIRECTORY CONTAINING BARCODE RNA AND DNA COUNTS 
#Read in all barcode counts and normalize by Reads per Million

filelist = list.files(path = '../../processed_data/expression_pipeline',
                      pattern = '^counts_*',
                      full.names = T)
for(i in filelist) {
    name <- gsub('counts_', '', basename(i))
    name <- gsub('.txt', '', name)
    x <- read.table(i, col.names=c(name, 'barcode'), header = F)
    x[[name]] <- 1000000*x[[name]]/sum(x[[name]])  #Normalizes by RPM
    assign(name,x)  
}

#combine reads for all barcodes 
rLP5_Endo2 <- full_join(rLP5_Endo2_DNA1, rLP5_Endo2_DNA2, by='barcode') %>%
  full_join(., rLP5_Endo2_RNA1, by='barcode') %>%
  full_join(., rLP5_Endo2_RNA2, by='barcode')

rm(rLP5_Endo2_DNA1, rLP5_Endo2_DNA2, rLP5_Endo2_RNA1, rLP5_Endo2_RNA2)
rm(x)

#Combine barcode counts with their promoter identity

barcode_stats_Endo2 <- read.table("../../processed_data/expression_pipeline/endo_mapping_barcode_statistics.txt", 
                                  header = T)
#Remove unmapped barcodes
mapped_barcodes <- barcode_stats_Endo2 %>% 
    filter(!is.na(most_common))
Compare_barcode_Reps <- left_join(mapped_barcodes, rLP5_Endo2 , by ='barcode') 
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0

temp <- filter(Compare_barcode_Reps, rLP5_Endo2_DNA1 > 0 | rLP5_Endo2_DNA2 > 0) #Remove barcodes with no DNA counts

#calculate expression of promoters
Endo2 <- temp %>% 
    group_by(most_common) %>% 
    mutate(num_barcodes = n()) %>%
    filter(rLP5_Endo2_DNA1 > 0 | rLP5_Endo2_DNA2 > 0) %>%
    mutate(num_barcodes_integrated = n()) %>%
    #Filter out promoters with fewer than 3 barcodes integrated
    filter(num_barcodes_integrated >= 4) %>% 
    mutate(DNA_sum = (sum(rLP5_Endo2_DNA2)+sum(rLP5_Endo2_DNA1)),
           RNA_exp_1 = sum(rLP5_Endo2_RNA1)/DNA_sum,
           RNA_exp_2 = sum(rLP5_Endo2_RNA2)/DNA_sum,
           RNA_exp_ave = ((RNA_exp_1+RNA_exp_2)/2)) %>% 
    ungroup() %>% 
    select(most_common, RNA_exp_1, RNA_exp_2, RNA_exp_ave, DNA_sum, num_barcodes, num_barcodes_integrated) %>% 
    distinct() 

#Convert sequence to variant
var_stats <- read.table("../../processed_data/expression_pipeline/endo_mapping_variant_statistics.txt",
                    header = T,
                    fill = T,
                    col.names = c('seq', 'name', 'num_barcodes', 
                                  'num_barcodes_unique', 'barcodes')) %>%
    select(name, seq) %>%
    mutate(name = gsub('>', '', name))
# remove any duplicates
var_stats_unique <- var_stats[!duplicated(var_stats$name), ]

Endo2 <- left_join(Endo2, var_stats_unique, by = c('most_common' = 'seq')) %>% 
    select('name', 'RNA_exp_1', 'RNA_exp_2', 'RNA_exp_ave', 'DNA_sum',
           'num_barcodes', 'num_barcodes_integrated') %>% 
    filter(!is.na(name))

# read in ref
ref <- read.table('../../ref/endo_lib_2016_controls_clean.txt', header = F, 
                  sep = '\t', col.names = c('name', 'seq')) %>% 
    mutate(name = gsub('>', '', name))

Endo2 <- left_join(Endo2, ref, by = 'name')

# remove primers
Endo2 <- Endo2 %>% 
    mutate(trimmed_seq = substring(seq, 25, 174))

# separate name into fields
Endo2 <- Endo2 %>% 
    separate(name, into = c('source', 'tss_pos', 'strand'), sep = ',', remove = F) %>% 
    mutate(tss_pos = as.numeric(tss_pos),
           strand = ifelse(is.na(strand), '+', strand))

# create accurate var start and end from TSS position
# parse negatives separately
neg <- subset(Endo2, grepl("neg_control", Endo2$name))

Endo2 <- Endo2 %>% 
    filter(!grepl('neg_control', name)) %>% 
    mutate(start = ifelse(strand == '+', tss_pos - 120, tss_pos - 30),
           end = ifelse(strand == '+', tss_pos + 30, tss_pos + 120))

Endo2 <- neg %>% 
    mutate(coord = gsub('neg_control_', '', name)) %>% 
    separate(coord, into = c('start', 'end'), sep = ':', convert = T) %>% 
    bind_rows(Endo2, .)

write.table(Endo2, "../../processed_data/expression_pipeline/rLP5_Endo2_expression.txt", quote = F, row.names = F)
