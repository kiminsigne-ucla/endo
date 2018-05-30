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


# run gkmSVM on TSS library sequences
# read input
active <- read.fasta('../processed_data/expression_pipeline/tss_positives.fasta')
# make sure names are unique
active_names <- paste0('seq', seq(1:length(active)))
# get train index
train_size = 0.75
train_index <- base::sample(seq(1:length(active)), size = length(active) * train_size, 
                            replace = F)
train <- active[train_index]
test <- active[-train_index]

write.fasta(train, names = active_names[train_index], nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/tss_positives_train.fasta')

write.fasta(test, names = active_names[-train_index], nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/tss_positives_test.fasta')

# read negative and split into train and test
negatives <- read.fasta('../processed_data/expression_pipeline/tss_negatives.fasta')
# get train index
neg_train_size = 0.75
neg_train_index <- base::sample(seq(1:length(negatives)), 
                                size = length(negatives) * neg_train_size, 
                                replace = F)
neg_train <- negatives[neg_train_index]
neg_test <- negatives[-neg_train_index]

write.fasta(neg_train, names = names(neg_train), nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/tss_negatives_train.fasta')

write.fasta(neg_test, names = names(neg_test), nbchar = 300,
            file.out = '../processed_data/gkmsvm_results/tss_negatives_test.fasta')


gkmsvm_kernel(posfile = "../processed_data/gkmsvm_results/tss_positives_train.fasta", 
              negfile = "../processed_data/gkmsvm_results/tss_negatives_train.fasta", 
              outfile = '../processed_data/gkmsvm_results/10mer_8ungapped_kernel', 
              L = 10, K = 8)

# perform SVM training with cross-validation
result <- gkmsvm_trainCV(kernelfn = '../processed_data/gkmsvm_results/10mer_8ungapped_kernel',
               posfn= "../processed_data/gkmsvm_results/tss_positives_train.fasta", 
               negfn = "../processed_data/gkmsvm_results/tss_negatives_train.fasta",
               svmfnprfx='../processed_data/gkmsvm_results/10mer_8ungapped_kernel', 
               outputCVpredfn='../processed_data/gkmsvm_results/10mer_8ungapped_cvpred',
               outputROCfn='../processed_data/gkmsvm_results/10mer_8ungapped_roc', 
               L = 10, K = 8,
               showPlots = T, outputPDFfn = '../processed_data/gkmsvm_results/10mer_8ungapped_curves')

ggsave('../processed_data/gkmsvm_results/10mer_8ungapped_ROC_PR_curves.png')


# # score all n-mers
# gkmsvm_classify(seqfile = '../processed_data/gkmsvm_results/all_10mers.fasta',
#                 svmfnprfx = '../processed_data/gkmsvm_results/promoter_svm_10mer_8ungapped_top10pct_300bp.out',
#                 outfile = '../processed_data/gkmsvm_results/all_10mers_classified_12mer_8ungapped_top10pct_300bp.txt',
#                 L = 10, K = 8)


