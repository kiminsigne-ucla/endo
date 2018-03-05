library(dplyr)
library(tidyr)
library(seqinr)
options(stringsAsFactors = F)

# Read in Data
Endo2 <- read.table("../processed_data/Endo2_integrated_alldata.txt", header = T, fill = T)

# read in ref
ref <- read.table('../ref/endo_lib_2016_controls_clean.txt', header = F, 
                  sep = '\t', col.names = c('name', 'seq')) %>% 
    mutate(name = gsub('>', '', name))

Endo2 <- left_join(Endo2, ref, by = 'name')

# identify 3standard deviations greater than median of negatives
neg_median <- median(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_average)
neg_sd <- sd(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_average)

# Subset all promoters that are greater than 3sd from the mean
positive_Endo2 <- filter(Endo2, RNA_exp_average > (neg_median+3*(neg_sd)))
negative_Endo2 <- filter(Endo2, RNA_exp_average < (neg_median+3*(neg_sd)))

positive_Endo2 %>% 
    separate(name, sep = ',', into = c('source', 'position', 'strand'), remove = F) %>% 
    mutate(chrom = 'U00096.2',
           start = as.numeric(position), 
           end = start + 1) %>% 
    select(chrom, start, end, name, RNA_exp_average, strand) %>% 
    mutate(strand = ifelse(grepl('neg', name), '+', strand)) %>% 
    filter(!is.na(start)) %>% 
    write.table('../processed_data/tss_positives.bed', sep = '\t', 
                col.names = F, row.names = F, quote = F)

negative_Endo2 %>% 
    separate(name, sep = ',', into = c('source', 'position', 'strand'), remove = F) %>% 
    mutate(chrom = 'U00096.2',
           start = ifelse(grepl('neg', name), var_start, as.numeric(position)), 
           end = ifelse(grepl('neg', name), var_end, start + 1)) %>% 
    select(chrom, start, end, name, RNA_exp_average, strand) %>% 
    mutate(strand = ifelse(grepl('neg', name), '+', strand)) %>% 
    filter(!is.na(start)) %>% 
    write.table('../processed_data/tss_negatives.bed', sep = '\t', 
                col.names = F, row.names = F, quote = F)

# write fasta
write.fasta(sequences = as.list(positive_Endo2$seq),
            names = positive_Endo2$name, nbchar = 200, as.string = T, open = 'w',
            file.out = '../processed_data/dragonn_output/tss_positives.fasta')

write.fasta(sequences = as.list(negative_Endo2$seq),
            names = negative_Endo2$name, nbchar = 200, as.string = T, open = 'w',
            file.out = '../processed_data/dragonn_output/tss_negatives.fasta')
