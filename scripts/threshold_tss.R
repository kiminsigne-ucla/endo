library(dplyr)
library(tidyr)
options(stringsAsFactors = F)

# Read in Data
Endo2 <- read.table("../processed_data/Endo2_integrated_alldata.txt", header = T, fill = T)

# identify 3standard deviations greater than median of negatives
neg_median <- median(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_average)
neg_sd <- sd(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_average)

# Subset all promoters that are greater than 3sd from the mean
positive_Endo2 <- filter(Endo2, RNA_exp_average > (neg_median+3*(neg_sd))) # > .7861676 == 2,048
negative_Endo2 <- filter(Endo2, RNA_exp_average < (neg_median+3*(neg_sd))) # < .7861676 == 15,215

positive_Endo2 %>% 
    separate(name, sep = ',', into = c('source', 'position', 'strand'), remove = F) %>% 
    mutate(chrom = 'U00096.2',
           start = as.numeric(position), 
           end = start + 1) %>% 
    select(chrom, start, end, name, RNA_exp_average, strand) %>% 
    filter(!is.na(start)) %>% 
    write.table('../processed_data/tss_positives.bed', sep = '\t', 
                col.names = F, row.names = F, quote = F)

negative_Endo2 %>% 
    separate(name, sep = ',', into = c('source', 'position', 'strand'), remove = F) %>% 
    mutate(chrom = 'U00096.2',
           start = ifelse(grepl('neg', name), var_start, as.numeric(position)), 
           end = ifelse(grepl('neg', name), var_end, start + 1)) %>% 
    select(chrom, start, end, name, RNA_exp_average, gene_strand) %>% 
    filter(!is.na(start)) %>% 
    write.table('../processed_data/tss_negatives.bed', sep = '\t', 
                col.names = F, row.names = F, quote = F)
