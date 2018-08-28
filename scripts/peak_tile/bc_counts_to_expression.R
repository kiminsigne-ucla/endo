#This script takes Barcode RNA and DNA counts, 
# merges them to their respective promoters, 
# and calculates the strength of each promoter


setwd("~/Documents/projects/ecoli_promoters/endo/scripts/peak_tile")
library("ggplot2")
library("dplyr")
library("tidyr")
require("cowplot")

options(stringsAsFactors = F)
options(scipen = 10000)

args = commandArgs(trailingOnly=TRUE)
args <- c('../../processed_data/peak_tile',
          '../../processed_data/peak_tile/peak_tile_mapping_barcode_statistics.txt',
          '../../processed_data/peak_tile/peak_tile_mapping_variant_statistics.txt',
          '../../ref/20180508_lb_peak_tile_lib.txt',
          '../../processed_data/peak_tile/peak_tile_expression.txt')
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
rLP5_peak_tile <- full_join(rLP5_PeakTiling_LB_DNA1_1_S5_R1_001, 
                               rLP5_PeakTiling_LB_DNA1_2_S6_R1_001, by='barcode') %>%
  full_join(., rLP5_PeakTiling_LB_DNA2_1_S7_R1_001, by='barcode') %>%
  full_join(., rLP5_PeakTiling_LB_DNA2_2_S8_R1_001, by='barcode') %>%
  full_join(., rLP5_PeakTiling_LB_RNA1_1_S1_R1_001, by='barcode') %>%
  full_join(., rLP5_PeakTiling_LB_RNA1_2_S2_R1_001, by='barcode') %>%
  full_join(., rLP5_PeakTiling_LB_RNA2_1_S3_R1_001, by='barcode') %>%
  full_join(., rLP5_PeakTiling_LB_RNA2_2_S4_R1_001, by='barcode')

names(rLP5_peak_tile) = c('DNA1_1', 'barcode', 'DNA1_2', 
                             'DNA2_1', 'DNA2_2', 
                             'RNA1_1', 'RNA1_2',
                             'RNA2_1', 'RNA2_2')
rLP5_peak_tile <- select(rLP5_peak_tile, barcode, DNA1_1, DNA1_2:RNA2_2)

rm(list = c('rLP5_PeakTiling_LB_DNA1_1_S5_R1_001',
            'rLP5_PeakTiling_LB_DNA1_2_S6_R1_001',
            'rLP5_PeakTiling_LB_DNA2_1_S7_R1_001',
            'rLP5_PeakTiling_LB_DNA2_2_S8_R1_001',
            'rLP5_PeakTiling_LB_RNA1_1_S1_R1_001',
            'rLP5_PeakTiling_LB_RNA1_2_S2_R1_001',
            'rLP5_PeakTiling_LB_RNA2_1_S3_R1_001',
            'rLP5_PeakTiling_LB_RNA2_2_S4_R1_001'))
rm(x)


#Combine barcode counts with their promoter identity
bc_stats <- read.table(bc_stats, header = T)

mapped_barcodes <- bc_stats[!is.na(bc_stats$most_common),] #Remove unmapped barcodes
Compare_barcode_Reps <- inner_join(mapped_barcodes, rLP5_peak_tile , by ='barcode') 
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0

nrow(Compare_barcode_Reps) / nrow(mapped_barcodes)

#Evaluate how many barcodes appear in integrated sequencing reads
nrow(filter(Compare_barcode_Reps, DNA1_1 > 0 | DNA2_1 > 0))

#Determine Expression as a function of the sum of all RNA counts for all barcodes 
#of a promoter divided by the sum of all DNA counts for all barcodes of a promoter
#Furthermore, filter out all promoters with fewer than 3 barcodes 

#temp <- filter(Compare_barcode_Reps, rLP5_peak_tile_LB_DNA1 > 0 | rLP5_peak_tile_LB_DNA2 > 0) #Remove barcodes with no DNA counts

peak_tile <- Compare_barcode_Reps %>% 
    group_by(most_common) %>% 
    mutate(num_barcodes_mapped = n()) %>%
    filter(DNA1_1 > 0, DNA1_2 > 0, DNA2_1 > 0, DNA2_2 > 0) %>%
    mutate(num_barcodes_integrated = n()) %>%
    filter(num_barcodes_integrated >= 3) %>% #Filter out promoters with fewer than 3 barcodes integrated
    mutate(RNA_exp_1_1 = sum(RNA1_1)/(sum(DNA1_1)),
           RNA_exp_1_2 = sum(RNA1_2)/(sum(DNA1_2)),
           RNA_exp_2_1 = sum(RNA2_1)/(sum(DNA2_1)),
           RNA_exp_2_2 = sum(RNA2_2)/(sum(DNA2_2)),
           RNA_exp_1 = ((RNA_exp_1_1+RNA_exp_1_2)/2),
           RNA_exp_2 = ((RNA_exp_2_1+RNA_exp_2_2)/2),
           RNA_exp_ave = ((RNA_exp_1+RNA_exp_2)/2),
           DNA_1 = (sum(DNA1_1)+sum(DNA1_2)), 
           DNA_2 = (sum(DNA2_1)+sum(DNA2_2)),  
           DNA_ave = ((DNA_1 + DNA_2)/2)) %>%
    ungroup() %>%
    select(most_common, RNA_exp_1_1, RNA_exp_1_2, RNA_exp_2_1, RNA_exp_2_2,
           RNA_exp_1, RNA_exp_2, RNA_exp_ave, DNA_1, DNA_2, DNA_ave, 
           num_barcodes_mapped, num_barcodes_integrated) %>% 
    distinct() 

var_stats <- read.table(var_stats, header = T, sep = '\t')
peak_tile <- var_stats %>%
    select(name, variant) %>% 
    inner_join(., peak_tile, by = c('variant' = 'most_common'))

write.table(peak_tile, output_name, quote = F, row.names = F)

percent_covered <- nrow(peak_tile) / nrow(var_stats)
# roughly 88% when requiring reads in all four DNA replicates


#distribution of the number of mapped and integrated barcodes
peak_tile %>% 
    select(mapped = num_barcodes_mapped, integrated = num_barcodes_integrated) %>% 
    gather('type', 'num_barcodes') %>% 
    ggplot(aes(num_barcodes)) + 
    # geom_density(aes(fill = type), alpha = 0.5) +
    geom_histogram(binwidth = 1, alpha = 0.5, position = 'identity', aes(fill = type)) + 
    labs(x = 'number of barcodes', fill = '')

#Plot expression correlations with negative controls and consensus promoters highlighted
test <- subset(peak_tile, grepl("neg_control", peak_tile$name))
test <- subset(test, !grepl("scramble", test$name))

corr <- summary(lm(RNA_exp_1 ~ RNA_exp_2, peak_tile))$r.squared

ggplot(peak_tile, aes(RNA_exp_1, RNA_exp_2)) + 
  geom_point(alpha = .4, aes(color = log2(DNA_ave))) +
  geom_point(data=test, aes(RNA_exp_1, RNA_exp_2), color = "firebrick1") +
  annotate("text", x =.1, y = 5, label = paste('R^2==', signif(corr, 3)), parse = T) + 
  annotation_logticks() + scale_x_log10(breaks = c(1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
  xlab('Expression (RNA/DNA) Replicate #1') + 
  ylab('Expression (RNA/DNA) Replicate #2') +
  viridis::scale_color_viridis() +
  ggtitle('Promoter Expression \n Between Biological Replicates') 

