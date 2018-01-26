source('https://bioconductor.org/biocLite.R')
biocLite('GenomicRanges')
biocLite('rtracklayer')
biocLite('BSgenome')
install.packages('ROCR')
install.packages('kernlab')
install.packages('seqinr')
install.packages('gkmSVM')

library('gkmSVM')
library('kernlab')
library('ROCR')
library('seqinr')
library('ggplot2')
library('dplyr')
require(cowplot)

options(stringsAsFactors = F)

# run gkmSVM on Basic-STARRseq calls
# 10mer/3 mismatches 
gkmsvm_kernel("./Enriched_150bp_.95qt_BasicStarrPeaks_train.fasta", 
              "./negative_150bp_.95qt_BasicStarrPeaks_train.fasta", 
              outfile = 'kernel.out')

# perform SVM training with cross-validation
gkmsvm_trainCV('kernel.out',
               # input files
               './Enriched_150bp_.95qt_BasicStarrPeaks_train.fasta',
               './negative_150bp_.95qt_BasicStarrPeaks_train.fasta', 
               svmfnprfx='promoter_svm.out', 
               outputCVpredfn='promoter_svm_cvpred.out',
               outputROCfn='promoter_svm_roc.out')
ggsave('ROC&P-R_Curves_10mer3mis.png')

gkmsvm_classify('../gkmSVM/nr10mers.fa',
                svmfnprfx='promoter_svm.out', 
                'promoter_svm_weights.out') #Score 10mers