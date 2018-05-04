#This script takes Barcode RNA and DNA counts, 
# merges them to their respective promoters, 
# and calculates the strength of each promoter

#Install required packages
# install.packages("dplyr")
# install.packages("ggplot2")

library("dplyr")

options(stringsAsFactors = F)

#SET TO WORKING DIRECTORY CONTAINING BARCODE RNA AND DNA COUNTS 
#Read in all barcode counts and normalize by Reads per Million

filelist = list.files(pattern = 'counts_*')
for(i in filelist) {
  x <- read.table(i, col.names=c(i, 'barcode'), header = F)
  x[[i]] <- 1000000*x[[i]]/sum(x[[i]])  #Normalizes by RPM
  assign(i,x)  
}
filelist

#combine reads for all barcodes 
rLP5_Endo2 <- full_join(rLP5_Endo2_DNA1.txt, rLP5_Endo2_DNA2.txt, by='barcode') %>%
  full_join(., rLP5_Endo2_RNA1.txt, by='barcode') %>%
  full_join(., rLP5_Endo2_RNA2.txt, by='barcode')

names(rLP5_Endo2) = sub(".txt","", names(rLP5_Endo2)) #rename all colummns that were named after text file
rm(list = c(filelist))
rm(x)



#Combine barcode counts with their promoter identity

barcode_stats_Endo2 <- read.table("./ref/barcode_statistics.txt", header = T)
mapped_barcodes <- barcode_stats_Endo2[!is.na(barcode_stats_Endo2$most_common),] #Remove unmapped barcodes
Compare_barcode_Reps <- left_join(mapped_barcodes, rLP5_Endo2 , by ='barcode') 
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0



temp <- filter(Compare_barcode_Reps, rLP5_Endo2_DNA1 > 0 | rLP5_Endo2_DNA2 > 0) #Remove barcodes with no DNA counts

#calculate expression of promoters
Endo2 <- temp %>% group_by(most_common) %>% 
  mutate(num_barcodes = n()) %>%
  filter(rLP5_Endo2_DNA1 > 0 | rLP5_Endo2_DNA2 > 0) %>%
  mutate(num_barcodes_integrated = n()) %>%
  filter(num_barcodes_integrated >= 4) %>% #Filter out promoters with fewer than 3 barcodes integrated
  mutate(RNA_exp_1 = sum(rLP5_Endo2_RNA1)/(sum(rLP5_Endo2_DNA1)+sum(rLP5_Endo2_DNA2)),
         RNA_exp_2 = sum(rLP5_Endo2_RNA2)/(sum(rLP5_Endo2_DNA2)+sum(rLP5_Endo2_DNA1)),
         RNA_exp_ave = ((RNA_exp_1+RNA_exp_2)/2),
         DNA_sum = (sum(rLP5_Endo2_DNA2)+sum(rLP5_Endo2_DNA1))) %>% 
  ungroup() %>% 
  select(most_common, RNA_exp_1, RNA_exp_2, RNA_exp_ave, DNA_sum, num_barcodes, num_barcodes_integrated) %>% 
  distinct() 

#Convert sequence to variant

Endo2 <- read.table("./ref/variant_statistics.txt",
                    header = T,
                    fill = T,
                    col.names = c('most_common', 'name', 'num_barcodes', 'num_barcodes_unique', 'barcodes')) %>%
  select(most_common, name) %>%
  inner_join(., Endo2, by = 'most_common') %>%
  distinct() %>%
  select('name', 'RNA_exp_1', 'RNA_exp_2', 'RNA_exp_ave', 'DNA_sum', 'num_barcodes', 'num_barcodes_integrated')

write.table(Endo2, "./rLP5_Endo2_expression.txt", quote = F, row.names = F)
Endo2 <- read.table("./rLP5_Endo2_expression.txt", header = T)

#Write out active and inactive promoters

neg_median <- median(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_ave)
neg_sd <- sd(subset(Endo2, grepl("neg_control", Endo2$name)) %>% .$RNA_exp_ave)

#read in coordinate files
coordinates <- read.table("./tss_coordinates.bed", header = F)
Endo2$name <- gsub('>', '', Endo2$name)

#Subset bed file for positive TSSs
tss_positives.bed <- filter(Endo2, RNA_exp_ave > (neg_median+2*(neg_sd))) %>%
                     semi_join(coordinates, . , by = c("V4" = "name")) %>%
		     write.table(tss_positives.bed, "./tss_positives.bed", col.names = F, row.names = F, quote = F, sep = '\t')

#Subset bed file for negative TSSs
tss_negatives.bed <- filter(Endo2, RNA_exp_ave < (neg_median+2*(neg_sd))) %>%
                     semi_join(coordinates, . , by = c("V4" = "name")) %>%
		     write.table(tss_negatives.bed, "./tss_negatives.bed", col.names = F, row.names = F, quote = F, sep = '\t')


