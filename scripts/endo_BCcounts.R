setwd("~/Documents/Kosuri Files/E. coli Genetics/Endo2.0/20170522_rLP5;Endo2/Datasets/")

library("ggplot2")
library("dplyr")
require(cowplot)

options(stringsAsFactors = F)

#Read in all Sequencing data for flp5;endo and normalize by Reads per Million

filelist = list.files(pattern = 'rLP5*')
for(i in filelist) {
  x <- read.table(i, col.names=c(i, 'barcode'), header = F)
  x[[i]] <- 1000000*x[[i]]/sum(x[[i]])  #Normalizes by RPM
  assign(i,x)  
}

#Combine Datasets
rLP5_Endo2 <- full_join(rLP5_Endo2_DNA1.txt, rLP5_Endo2_DNA2.txt, by='barcode') %>%
  full_join(., rLP5_Endo2_RNA1.txt, by='barcode') %>%
  full_join(., rLP5_Endo2_RNA2.txt, by='barcode') 
  

names(rLP5_Endo2) = sub(".txt","", names(rLP5_Endo2)) #rename all colummns that were named after text file temp <- Loop_Data[c(-1,-5)]
rm(list = c(filelist))
rm(x)

#identify barcodes with mapped variants
barcode_stats <- read.table("./barcode_statistics.txt", header = T)

Compare_barcode_Reps <- left_join(barcode_stats, rLP5_Endo2 , by ='barcode')
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0

#Count how many mapped barcodes appear *optional*
nrow(filter(Compare_barcode_Reps, rLP5_Endo2_DNA1 > 0 & rLP5_Endo2_DNA2 > 0)) #1,072,014/2,011,019 (53.3% of integrants)


#Calculate promoter Expression, filtering out barcodes without DNA reads and filtering out variants with less than 6 barcodes per variant

Promoter_Expression <- Compare_barcode_Reps %>% group_by(name) %>% 
  mutate(num_mapped_barcodes = n()) %>%
  filter(rLP5_Endo2_DNA1 > 0 & rLP5_Endo2_DNA2 > 0) %>%
  mutate(num_integrated_barcodes = n()) %>%
  filter(num_integrated_barcodes >= 6) %>%
  mutate(RNA_exp_1 = (sum(rLP5_Endo2_RNA1)/sum(rLP5_Endo2_DNA1)),
         RNA_exp_2 = (sum(rLP5_Endo2_RNA2)/sum(rLP5_Endo2_DNA2)),
         DNA_sum_2 = sum(rLP5_Endo2_DNA2), 
         DNA_sum_1 = sum(rLP5_Endo2_DNA1), 
         RNA_exp_average = (((RNA_exp_1) + (RNA_exp_2))/2),
         DNA_sum_total = (DNA_sum_1 + DNA_sum_2)) %>% 
  ungroup() %>% 
  select(name, RNA_exp_1, RNA_exp_2,RNA_exp_average, DNA_sum_1, DNA_sum_2,DNA_sum_total, num_mapped_barcodes, num_integrated_barcodes) %>% 
  distinct() 

write.table(Promoter_Expression, "./Endo2_Promoter_expression.txt", row.names = F, quote = F)

#Plot correlation between technical replicates, highlighting 
test <- subset(Promoter_Expression, grepl("neg_control", Promoter_Expression$name))

corr <- summary(lm(RNA_exp_1 ~ RNA_exp_2, Promoter_Expression))$r.squared

ggplot(Promoter_Expression, aes(RNA_exp_1, RNA_exp_2)) + geom_point() +
  geom_smooth(method=lm) + annotate("text", x =.5, y = 30, label = paste('r^2==', signif(corr, 3)), parse = T) +
  xlim(0,100) +  scale_x_log10() + scale_y_log10() +
  xlab('LB Expression Rep1') + ylab('LB Expression Rep2') + ggtitle('Comparing Promoter Expression Between Technical Replicates') +
  geom_point(data=test, aes(RNA_exp_1, RNA_exp_2), color = "firebrick1")
ggsave('../Figures/CorrelationBetweenTechnicalReplicates.pdf')


#Plot histogram of integrated barcodes per variant
ggplot(Promoter_Expression, aes(num_mapped_barcodes)) + geom_density(alpha=.2, fill = "#3B9AB2", color = "#3B9AB2") +
  geom_density(aes(num_integrated_barcodes), alpha=.2, fill = "#F21A00", color = "#F21A00") +
  scale_x_log10() + annotation_logticks(sides = 'b') +
  annotate("text", x=13, y=1.35, label = paste(' Integrated Median =', median(Promoter_Expression$num_integrated_barcodes), 'barcodes')) +
  annotate("text", x=13, y=1.43, label = paste('Mapped Median =', median(Promoter_Expression$num_mapped_barcodes), 'barcodes')) +
  labs(x='Number of Mapped Barcodes', y='Density', title = 'Distribution of the Number of Barcodes Per Variant')
ggsave('../Figures/Distribution of Barcodes per Variant.pdf')

#add coordinates to promoter

#var_annot <- read.table("./var_annot.txt", header = T, fill = T)
#var_coords <- var_annot[c(1,as.numeric(25),as.numeric(26),32)]
#var_coords$name <- gsub(">", "", var_coords$name)
#write.table(var_coords, "./endo2VariantCoordinates.txt", row.names = F, quote = F)
var_coords <- read.table("./endo2VariantCoordinates.txt", header = T, fill = T)

Endo2_integrated_alldata <- left_join(Promoter_Expression, var_coords, by = 'name')
Endo2_integrated_alldata$var_start <- as.numeric(as.character(Endo2_integrated_alldata$var_start))
Endo2_integrated_alldata$var_end <- as.numeric(as.character(Endo2_integrated_alldata$var_end))

write.table(Endo2_integrated_alldata, "./Endo2_integrated_alldata.txt", row.names = F, quote = F)












       