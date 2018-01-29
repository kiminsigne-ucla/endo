# run peak calling in R
Rscript run_basic_starrseq.R 300 0.90 top10pct_300bp

# # format results
# ./format_basic_starrseq.sh ../processed_data/basic_starrseq_results/plus_peaks_top5pct_100bp.bed \
# ../processed_data/basic_starrseq_results/minus_peaks_top5pct_100bp.bed \
# ../processed_data/basic_starrseq_results/plus_minus_peaks_top5pct_100bp.fasta \
# ../processed_data/basic_starrseq_results/negative_peak_controls_top5pct.fasta
