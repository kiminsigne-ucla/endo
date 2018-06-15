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

# read input
peaks <- read.fasta('../processed_data/frag_peak_calling/U00096.2_plus_minus_called_peaks_threshold1.1_merge40_min60.fasta')

# get train index
train_size = 0.75
train_index <- base::sample(seq(1:length(peaks)), 
                            size = length(peaks) * train_size, 
                            replace = F)
train <- peaks[train_index]
test <- peaks[-train_index]

write.fasta(train, names = names(train), nbchar = 5000,
            file.out = '../processed_data/gkmsvm_results/peaks/train.fasta')

write.fasta(test, names = names(test), nbchar = 5000,
            file.out = '../processed_data/gkmsvm_results/peaks/test.fasta')

# read negative and split into train and test
negatives <- read.fasta('../processed_data/frag_peak_calling/U00096.2_negative_generated_peaks_threshold1.1_merge40_min60.fasta')
names(negatives) <- paste0('neg_', names(negatives))
# get train index
neg_train_size = 0.75
neg_train_index <- base::sample(seq(1:length(negatives)), 
                                size = length(negatives) * neg_train_size, 
                                replace = F)
neg_train <- negatives[neg_train_index]
neg_test <- negatives[-neg_train_index]

write.fasta(neg_train, names = names(neg_train), nbchar = 5000,
            file.out = '../processed_data/gkmsvm_results/peaks/negative_train.fasta')

write.fasta(neg_test, names = names(neg_test), nbchar = 5000,
            file.out = '../processed_data/gkmsvm_results/peaks/negative_test.fasta')


gkmsvm_kernel(posfile = "../processed_data/gkmsvm_results/peaks/train.fasta", 
              negfile = "../processed_data/gkmsvm_results/peaks/negative_train.fasta", 
              outfile = '../processed_data/gkmsvm_results/peaks/10mer_8ungapped_kernel', 
              L = 10, K = 8)

# perform SVM training with cross-validation
gkmsvm_trainCV(kernelfn = '../processed_data/gkmsvm_results/peaks/10mer_8ungapped_kernel',
               posfn= "../processed_data/gkmsvm_results/peaks/train.fasta", 
               negfn = "../processed_data/gkmsvm_results/peaks/negative_train.fasta",
               svmfnprfx='../processed_data/gkmsvm_results/peaks/frag_peak_svm', 
               outputCVpredfn='../processed_data/gkmsvm_results/peaks/frag_peak_cvpred_10mer_8ungapped',
               outputROCfn='../processed_data/gkmsvm_results/peaks/frag_peak_roc_10mer_8ungapped', 
               L = 10, K = 8)

ggsave('../processed_data/gkmsvm_results/peaks/frag_peak_10mer_8ungapped_ROC_PR_curves.png')

# classify test
gkmsvm_classify(seqfile = '../processed_data/gkmsvm_results/peaks/test.fasta',
                svmfnprfx = '../processed_data/gkmsvm_results/peaks/frag_peak_svm',
                outfile = '../processed_data/gkmsvm_results/peaks/frag_peak_svm_test.txt',
                L = 10, K = 8)

gkmsvm_classify(seqfile = '../processed_data/gkmsvm_results/peaks/negative_test.fasta',
                svmfnprfx = '../processed_data/gkmsvm_results/peaks/frag_peak_svm',
                outfile = '../processed_data/gkmsvm_results/peaks/frag_peak_svm_negative_test.txt',
                L = 10, K = 8)

results <- read.table('../processed_data/gkmsvm_results/peaks/frag_peak_svm_test.txt',
                      header = F)
colnames(results) <- c('name', 'score')
results <- results %>%
    mutate(type = 'positive')

neg_results <- read.table('../processed_data/gkmsvm_results/peaks/frag_peak_svm_negative_test.txt',
                          header = F)
colnames(neg_results) <- c('name', 'score')
results <- neg_results %>%
    mutate(type = 'negative') %>%
    bind_rows(results)

library(pROC)
roc_obj <- roc(results$type, results$score)
auroc <- auc(roc_obj)
roc_df <- data.frame(
    tpr = rev(roc_obj$sensitivities),
    fpr = rev(1 - roc_obj$specificities))

ggplot(roc_df, aes(fpr, tpr)) + geom_line() +
    annotate('text', x = 0.75, y = 0.10, label = paste0("AUC = ", signif(auroc, 3))) +
    labs(x = '1 - specificity', y = 'sensitivity', 
         title = 'gkm-SVM test set (n = 1740)')

ggplot(results, aes(score)) + geom_histogram(binwidth = 0.10, aes(fill = type))

ggplot(results, aes(score)) + geom_density(aes(color = type))

ggplot(results, aes(score)) + stat_ecdf(aes(color = type)) +
    labs(x = 'SVM score', y = '', title = 'CDF of SVM scores')
