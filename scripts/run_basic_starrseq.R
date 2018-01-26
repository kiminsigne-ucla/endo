#https://bioconductor.org/packages/release/bioc/html/BasicSTARRseq.html

#Use starrseq to call peaks from data
source("https://bioconductor.org/biocLite.R")
# biocLite("BasicSTARRseq")
library(BasicSTARRseq)
library(dplyr)
#install.packages(matricstats))


run_basic_starrseq <- function(rna_file, dna_file, peak_width = 100, quantile = 0.90) {
    
    starrseqGRanges <- granges(readGAlignments(rna_file)) 
    inputGRanges <- granges(readGAlignments(dna_file)) 
    data <- STARRseqData(sample = starrseqGRanges, control = inputGRanges, pairedEnd = F) 
    data
    
    peaks <- getPeaks(data, deduplicate = F, peakWidth = peak_width, minQuantile = quantile)
    head(peaks)
    
    peaks_df<- as.data.frame(peaks)
    
    return(peaks_df)
}

write_results <- function(peaks_df, new_strand, output_file) {
    # format: chr, start, end, name, score, strand
    # name of each peak will be format start,end,strand
    # score will be format pVal, enrichment
    peaks_df %>% 
        mutate(strand = new_strand,
               name = paste(start, end, strand, sep = ','),
               score = paste(pVal, enrichment, sep = ',')) %>% 
        select(seqnames, start, end, name, score, strand) %>% 
        write.table(output_file, quote = F, row.names = F, col.names = F, sep = '\t')
}

peaks_minus <- run_basic_starrseq("../rawdata/genome_frag/Minus_RNA.bam", 
                                  "../rawdata/genome_frag/Minus_DNA.bam",
                                  peak_width = 100, quantile = 0.95)
write_results(peaks_minus, 
              new_strand = '-',
              "../processed_data/basic_starrseq_results/minus_peaks_top5pct_100bp.bed")

peaks_plus <- run_basic_starrseq("../rawdata/genome_frag/Plus_RNA.bam",
                                 "../rawdata/genome_frag/Plus_DNA.bam",
                                 peak_width = 100, quantile = 0.95) 
write_results(peaks_plus,
              new_strand = '+',
              "../processed_data/basic_starrseq_results/plus_peaks_top5pct_100bp.bed")

