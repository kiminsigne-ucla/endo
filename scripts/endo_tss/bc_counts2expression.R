#This script takes Barcode RNA and DNA counts, 
# merges them to their respective promoters, 
# and calculates the strength of each promoter

#Install required packages
# install.packages("dplyr")
# install.packages("ggplot2")

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
options(stringsAsFactors = F)

# expect five arguments in this order: relative path to barcode count files,
# barcode statistics file, variant statistics file, library reference file,
# output name
args = commandArgs(trailingOnly=TRUE)
args <- c('../../processed_data/expression_pipeline',
          '../../processed_data/expression_pipeline/endo_mapping_imperfect_barcode_statistics.txt',
          '../../processed_data/expression_pipeline/endo_mapping_imperfect_variant_statistics.txt',
          '../../processed_data/expression_pipeline/endo_lib_imperfects.txt',
          '../../processed_data/expression_pipeline/endo_imperfect_expression.txt')
count_folder <- args[1]
bc_stats <- args[2]
var_stats <- args[3]
lib_file <- args[4]
output_name <- args[5]

#SET TO WORKING DIRECTORY CONTAINING BARCODE RNA AND DNA COUNTS 
#Read in all barcode counts and normalize by Reads per Million

filelist = list.files(path = count_folder,
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

barcode_stats_Endo2 <- read.table(bc_stats, header = T)
#Remove unmapped barcodes
mapped_barcodes <- barcode_stats_Endo2 %>% 
    filter(!is.na(most_common))
Compare_barcode_Reps <- left_join(mapped_barcodes, rLP5_Endo2 , by ='barcode') 
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0

#Remove barcodes with no DNA counts
temp <- filter(Compare_barcode_Reps, rLP5_Endo2_DNA1 > 0 | rLP5_Endo2_DNA2 > 0) 

#calculate expression of promoters
barcode_cutoff <- 1
Endo2 <- temp %>% 
    filter(rLP5_Endo2_DNA1 > 0 | rLP5_Endo2_DNA2 > 0) %>%
    group_by(most_common) %>% 
    mutate(num_barcodes = n(),
           num_barcodes_integrated = n()) %>%
    #Filter out promoters with fewer than 3 barcodes integrated
    filter(num_barcodes_integrated >= barcode_cutoff) %>% 
    mutate(DNA_sum = (sum(rLP5_Endo2_DNA2)+sum(rLP5_Endo2_DNA1)),
           RNA_exp_1 = sum(rLP5_Endo2_RNA1)/DNA_sum,
           RNA_exp_2 = sum(rLP5_Endo2_RNA2)/DNA_sum,
           RNA_exp_ave = ((RNA_exp_1+RNA_exp_2)/2)) %>% 
    ungroup() %>% 
    select(most_common, RNA_exp_1, RNA_exp_2, RNA_exp_ave, DNA_sum, num_barcodes, num_barcodes_integrated) %>% 
    distinct()

#Convert sequence to variant
var_stats <- read.table(var_stats,
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
           'num_barcodes', 'num_barcodes_integrated', 'most_common') %>% 
    filter(!is.na(name))

# read in ref
ref <- read.table(lib_file, header = F, 
                  sep = '\t', col.names = c('name', 'seq')) %>% 
    mutate(name = gsub('>', '', name))

Endo2 <- left_join(Endo2, ref, by = 'name')


if(!grepl('imperfect', output_name)){
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
    
    # get original (unflipped) sequence based on coordinates
    Endo2 %>% 
        mutate(chrom = 'U00096.2',
               score = '.') %>% 
        select(chrom, start, end, name, score, strand) %>%
        filter(start > 0) %>% 
        write.table(file = '../../processed_data/expression_pipeline/endo_bed.txt',
                    sep = '\t', quote = F, col.names = F, row.names = F)
    
    system(paste(' bedtools getfasta',
                 '-fi ../../ref/U00096.2.fasta',
                 '-bed ../../processed_data/expression_pipeline/endo_bed.txt',
                 '-fo ../../ref/endo_lib_original_seq.txt',
                 '-name -s -tab'))
    
    original <- read.table('../../ref/endo_lib_original_seq.txt', header = F,
                           col.names = c('name', 'original_seq'))
    
    Endo2 <- left_join(Endo2, original, by = 'name') %>% 
        select(-seq)
    
    # write simple file
    Endo2 %>% 
        select(original_seq, RNA_exp_ave) %>% 
        filter(!is.na(original_seq)) %>% 
        write.table(file = '../../processed_data/expression_pipeline/tss_all.txt',
                    sep = '\t',
                    quote = F, row.names = F, col.names = F)
}

write.table(Endo2, file = output_name, quote = F, row.names = F)

# Endo2 <- neg %>% 
#     mutate(coord = gsub('neg_control_', '', name)) %>% 
#     separate(coord, into = c('start', 'end'), sep = ':', convert = T) %>% 
#     bind_rows(Endo2, .)
# 
# Endo2 <- Endo2 %>% 
#     mutate(exp_sum = RNA_exp_1 + RNA_exp_2)
# 
# write.table(Endo2, output_name, quote = F, row.names = F)
neg <- subset(Endo2, grepl("neg_control", Endo2$name))
neg_median <- median(neg$RNA_exp_ave)
neg_sd <- sd(neg$RNA_exp_ave)

corr <- cor(Endo2$RNA_exp_1, Endo2$RNA_exp_2)   
Endo2 %>% 
    mutate(category = ifelse(grepl('neg_control', name), 'negative', 'tss')) %>% 
    ggplot(aes(RNA_exp_1, RNA_exp_2)) + geom_point(aes(color = category)) +
        scale_x_log10() + scale_y_log10() + annotation_logticks(sides = 'bl') +
        geom_vline(xintercept = neg_median + 3*neg_sd, linetype = 'dashed') +
        scale_color_manual(values = c('red', 'black')) + 
        labs(x = 'replicate 1', y = 'replicate 2', color = '') +
        annotate('text', x = 10, y = 0.1, label = paste0('r = ', signif(corr, 3))) +
        theme(legend.position = 'none')

ggplot(Endo2, aes(num_barcodes_integrated)) + geom_histogram(binwidth = 1) +
    labs(x = 'number of integrated barcodes') +
    geom_vline(xintercept = median(Endo2$num_barcodes_integrated))
