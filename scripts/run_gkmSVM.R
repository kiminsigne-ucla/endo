# source('https://bioconductor.org/biocLite.R')
# biocLite('GenomicRanges')
# biocLite('rtracklayer')
# biocLite('BSgenome')
# install.packages('ROCR')
# install.packages('kernlab')
# install.packages('seqinr')
# install.packages('gkmSVM')

library('gkmSVM')
library('kernlab')
library('ROCR')
library('seqinr')
library('ggplot2')
library('dplyr')
require(cowplot)

options(stringsAsFactors = F)
set.seed(123)


# run gkmSVM on BasicSTARRseq calls
# read input
peaks <- read.fasta('../processed_data/basic_starrseq_results/plus_minus_peaks_top10pct_300bp.fasta')
# make sure names are unique
peak_names <- paste0('seq', seq(1:length(peaks)))
# get train index
train_size = 0.75
train_index <- base::sample(seq(1:length(peaks)), size = length(peaks) * train_size, 
                      replace = F)
train <- peaks[train_index]
test <- peaks[-train_index]

write.fasta(train, names = peak_names[train_index], nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/plus_minus_peaks_top10pct_300bp_train.fasta')

write.fasta(test, names = peak_names[-train_index], nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/plus_minus_peaks_top10pct_300bp_test.fasta')

# read negative and split into train and test
negatives <- read.fasta('../processed_data/basic_starrseq_results/negative_peaks_controls_top10pct_300bp.fasta')
# get train index
neg_train_size = 0.75
neg_train_index <- base::sample(seq(1:length(negatives)), 
                          size = length(negatives) * neg_train_size, 
                          replace = F)
neg_train <- negatives[neg_train_index]
neg_test <- negatives[-neg_train_index]

write.fasta(neg_train, names = names(neg_train), nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/negative_peaks_controls_top10pct_300bp_train.fasta')

write.fasta(neg_test, names = names(neg_test), nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/negative_peaks_controls_top10pct_300bp_test.fasta')

# 12mer, 8ungapped nucleotides
gkmsvm_kernel(posfile = "../processed_data/gkmsvm_results/plus_minus_peaks_top10pct_300bp_train.fasta", 
              negfile = "../processed_data/gkmsvm_results/negative_peaks_controls_top10pct_300bp_train.fasta", 
              outfile = '../processed_data/gkmsvm_results/12mer_8ungapped_kernel_300bp_top10pct.out', 
              L = 12, K = 8)

# perform SVM training with cross-validation
gkmsvm_trainCV(kernelfn = '../processed_data/gkmsvm_results/12mer_8ungapped_kernel_300bp_top10pct.out',
               posfn= "../processed_data/gkmsvm_results/plus_minus_peaks_top10pct_300bp_train.fasta", 
               negfn = "../processed_data/gkmsvm_results/negative_peaks_controls_top10pct_300bp_train.fasta",
               svmfnprfx='../processed_data/gkmsvm_results/promoter_svm_12mer_8ungapped_top10pct_300bp.out', 
               outputCVpredfn='../processed_data/gkmsvm_results/promoter_svm_cvpred_12mer_8ungapped_top10pct_300bp.out',
               outputROCfn='../processed_data/gkmsvm_results/promoter_svm_roc_12mer_8ungapped_top10pct_300bp.out', 
               L = 12, K = 8)

ggsave('../processed_data/gkmsvm_results/300bp_top10pct_12mer_8ungapped_ROC_PR_curves.png')


# score all 12-mers
gkmsvm_classify(seqfile = '../processed_data/gkmsvm_results/all_12mers.fasta',
                svmfnprfx = '../processed_data/gkmsvm_results/promoter_svm_12mer_8ungapped_top10pct_300bp.out',
                outfile = '../processed_data/gkmsvm_results/all_12mers_classified_12mer_8ungapped_top10pct_300bp.txt',
                L = 12, K = 8)

