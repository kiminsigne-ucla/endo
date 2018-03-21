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

# remove primers
Endo2 <- Endo2 %>% 
    mutate(trimmed_seq = substring(seq, 25, 174))

# separate name into fields
Endo2 <- Endo2 %>% 
    separate(name, into = c('source', 'tss_pos', 'strand'), sep = ',', remove = F)

Endo2 <- Endo2 %>% 
    mutate(tss_pos = as.numeric(tss_pos),
           strand = ifelse(is.na(strand), '+', strand))

# identify 3standard deviations greater than median of negatives
neg_median <- median(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_average)
neg_sd <- sd(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_average)

# create accurate var start and end from TSS position
Endo2 <- Endo2 %>% 
    mutate(start = ifelse(strand == '+', tss_pos - 120, tss_pos - 30),
           end = ifelse(strand == '+', tss_pos + 30, tss_pos + 120)) %>% 
    mutate(start = ifelse(grepl('neg', name), var_start, start),
           end = ifelse(grepl('neg', name), var_end, end))

# Subset all promoters that are greater than 3sd from the mean
positive_Endo2 <- filter(Endo2, RNA_exp_average > (neg_median+3*(neg_sd)))
negative_Endo2 <- filter(Endo2, RNA_exp_average < (neg_median+3*(neg_sd)))

false_negatives <- read.table('../processed_data/tss_false_negative.bed', sep='\t',
                              col.names = c('chrom', 'start', 'end', 'name'))

pos_with_false_neg <- bind_rows(positive_Endo2, 
                                semi_join(Endo2, false_negatives, by = 'name'))

neg_without_false_neg <- anti_join(negative_Endo2, false_negatives, by = 'name')

# positive_Endo2 %>% 
#     separate(name, sep = ',', into = c('source', 'position', 'strand'), remove = F) %>% 
#     mutate(chrom = 'U00096.2',
#            start = as.numeric(position), 
#            end = start + 1) %>% 
#     select(chrom, start, end, name, RNA_exp_average, strand) %>% 
#     mutate(strand = ifelse(grepl('neg', name), '+', strand)) %>% 
#     filter(!is.na(start)) %>% 
#     write.table('../processed_data/tss_positives.bed', sep = '\t', 
#                 col.names = F, row.names = F, quote = F)

positive_Endo2 %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_average, strand) %>%
    filter(start > 0) %>%
    write.table('../processed_data/tss_positives.bed', sep = '\t',
                col.names = F, row.names = F, quote = F)

pos_with_false_neg %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_average, strand) %>%
    filter(start > 0) %>%
    write.table('../processed_data/tss_positives_with_false_neg.bed', sep = '\t',
                col.names = F, row.names = F, quote = F)

# negative_Endo2 %>% 
#     separate(name, sep = ',', into = c('source', 'position', 'strand'), remove = F) %>% 
#     mutate(chrom = 'U00096.2',
#            start = ifelse(grepl('neg', name), var_start, as.numeric(position)), 
#            end = ifelse(grepl('neg', name), var_end, start + 1)) %>% 
#     select(chrom, start, end, name, RNA_exp_average, strand) %>% 
#     mutate(strand = ifelse(grepl('neg', name), '+', strand)) %>% 
#     filter(!is.na(start)) %>% 
#     write.table('../processed_data/tss_negatives.bed', sep = '\t', 
#                 col.names = F, row.names = F, quote = F)

negative_Endo2 %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_average, strand) %>%
    filter(start > 0) %>%
    write.table('../processed_data/tss_negatives.bed', sep = '\t',
                col.names = F, row.names = F, quote = F)

neg_without_false_neg %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_average, strand) %>%
    filter(start > 0) %>%
    write.table('../processed_data/tss_negatives_without_false_neg.bed', sep = '\t',
                col.names = F, row.names = F, quote = F)

# # write fasta
# write.fasta(sequences = as.list(positive_Endo2$trimmed_seq),
#             names = positive_Endo2$name, nbchar = 200, as.string = T, open = 'w',
#             file.out = '../processed_data/tss_positives.fasta')
# 
# write.fasta(sequences = as.list(negative_Endo2$trimmed_seq),
#             names = negative_Endo2$name, nbchar = 200, as.string = T, open = 'w',
#             file.out = '../processed_data/tss_negatives.fasta')

Endo2 %>% 
    left_join(seqs, by = 'name') %>% 
    select(original_seq, RNA_exp_average) %>% 
    write.table(file = '../processed_data/tss_all.txt', sep = '\t', 
                col.names = F, row.names = F, quote = F)

