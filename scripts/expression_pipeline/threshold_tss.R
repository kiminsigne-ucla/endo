Endo2 <- read.table("rLP5_Endo2_expression.txt", header = T)

# Write out active and inactive promoters
positive_Endo2 %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_ave, strand) %>%
    filter(start > 0) %>%
    write.table('../../processed_data/expression_pipeline/tss_positives.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

negative_Endo2 %>%
    mutate(chrom = 'U00096.2') %>%
    select(chrom, start, end, name, RNA_exp_ave, strand) %>%
    filter(start > 0) %>%
    write.table('../../processed_data/expression_pipeline/tss_negatives.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)