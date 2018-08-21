library("dplyr")
options(stringsAsFactors = F)

#Read in all Sequencing data for rLP5-Frag

filelist = list.files(path = '../../processed_data/frag_peak_calling',
                      pattern = 'rLP5_frag*',
                      full.names = T)
for(i in filelist) {
    name <- gsub('.txt', '',basename(i))
    x <- read.table(i, col.names=c(name, 'barcode'), header = F, fill=T)
    x[[name]] <- 1000000*x[[name]]/sum(x[[name]])  #Normalizes by RPM
    assign(name,x)
}

rLP5_frag <- full_join(rLP5_frag_DNA1_1, rLP5_frag_DNA1_2, by='barcode') %>%
  full_join(rLP5_frag_DNA2_1, by='barcode') %>%
  full_join(rLP5_frag_DNA2_2, by='barcode') %>%
  full_join(rLP5_frag_RNA1_1, by='barcode') %>%
  full_join(rLP5_frag_RNA1_2, by='barcode') %>%
  full_join(rLP5_frag_RNA2_1, by='barcode') %>%
  full_join(rLP5_frag_RNA2_2, by='barcode') 

rm(rLP5_frag_DNA1_1, rLP5_frag_DNA1_2, rLP5_frag_DNA2_1, rLP5_frag_DNA2_2,
   rLP5_frag_RNA1_1, rLP5_frag_RNA1_2, rLP5_frag_RNA2_1, rLP5_frag_RNA2_2)
rm(x)

rLP5_frag[is.na(rLP5_frag)] <- 0

mapped_fragstats <- read.table("../../processed_data/frag_peak_calling/frag_stats.txt", 
                               header = F, fill = T, 
                               col.names = c('barcode', 'start', 'end', 'strand', 'fragment'))
mapped_fragstats <- mapped_fragstats[!duplicated(mapped_fragstats$barcode),]

fragstats <- left_join(mapped_fragstats, rLP5_frag, by ='barcode')

Frag_LB <- fragstats %>%
  group_by(fragment) %>% 
  mutate(num_mapped_barcodes = n()) %>%
  filter(rLP5_frag_DNA1_1 > 0 | rLP5_frag_DNA1_2 > 0 | rLP5_frag_DNA2_1 > 0 | rLP5_frag_DNA2_2 > 0) %>%
  mutate(num_integrated_barcodes = n()) %>%
  filter(num_integrated_barcodes >= 1) %>%
  mutate(RNA_exp_1 = (sum(rLP5_frag_RNA1_1)+sum(rLP5_frag_RNA1_2))/(sum(rLP5_frag_DNA1_1)+sum(rLP5_frag_DNA1_2)),
         RNA_exp_2 = (sum(rLP5_frag_RNA2_1)+sum(rLP5_frag_RNA2_2))/(sum(rLP5_frag_DNA2_1)+sum(rLP5_frag_DNA2_2)),
         RNA_exp_ave = ((RNA_exp_1 + RNA_exp_2)/2),
         DNA_sum_1 = sum(rLP5_frag_DNA1_1)+sum(rLP5_frag_DNA1_2),
         DNA_sum_2 = sum(rLP5_frag_DNA1_1)+sum(rLP5_frag_DNA1_2),
         DNA_ave = ((DNA_sum_1 + DNA_sum_2)/2)) %>% 
  ungroup() %>% 
  filter(RNA_exp_1 >= .1 & RNA_exp_2 >= .1, RNA_exp_1 < 10000 & RNA_exp_2 < 10000, DNA_ave>1) %>% #Filter out fragments with low expression of 
  mutate(variation = abs(log2(RNA_exp_2/RNA_exp_1))) %>% 
  filter(variation < 2.32) %>% #Remove fragments that have a 5x difference between biological replicates
  select(fragment, RNA_exp_1, RNA_exp_2, RNA_exp_ave, DNA_sum_1, DNA_sum_2, DNA_ave, num_mapped_barcodes, num_integrated_barcodes, start, end, strand, variation) %>% 
  distinct() 

write.table(Frag_LB, "../../processed_data/frag_peak_calling/U00096.2_frag-rLP5_LB_expression.txt", quote = F, row.names = F)

Frag_LB %>% 
    select(fragment, RNA_exp_ave) %>% 
    write.table('../../processed_data/frag_peak_calling/frag_expression.txt', sep = '\t',
                quote = F, row.names = F, col.names = F)

Frag_LB <- read.table('../../processed_data/frag_peak_calling/U00096.2_frag-rLP5_LB_expression.txt',
                      header = T) %>% 
    mutate(frag_length = end - start)

ggplot(Frag_LB, aes(frag_length)) + geom_histogram(binwidth = 1) +
    geom_vline(xintercept = median(Frag_LB$frag_length), linetype = 'dashed') +
    labs(x = 'fragment length')

corr <- cor(Frag_LB$RNA_exp_1, Frag_LB$RNA_exp_2)
ggplot(Frag_LB, aes(RNA_exp_1, RNA_exp_2)) + geom_point(alpha = 0.50) +
    scale_x_log10() + scale_y_log10() + annotation_logticks(sides = 'bl') +
    labs(x = 'replicate 1', y = 'replicate 2') +
    annotate('text', x = 100, y = 0.5, label = paste0('r = ', signif(corr, 3)),
             size = 6)

peaks <- read.table('../../processed_data/frag_peak_calling/U00096.2_plus_minus_called_peaks_threshold1.1_merge40_min60.bed',
                    col.names = c('chrom', 'start', 'end', 'name', 'total_exp', 'strand'))
peaks <- peaks %>% 
    mutate(length = end - start)

ggplot(peaks, aes(length)) + geom_histogram(binwidth = 10) +
    geom_vline(xintercept = median(peaks$length), linetype = 'dashed') +
    labs(x = 'peak length')
