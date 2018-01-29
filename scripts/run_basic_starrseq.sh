# run peak calling in R
Rscript run_basic_starrseq.R 

# format results
./format_basic_starrseq.sh ../processed_data/basic_starrseq_results/plus_peaks_top10pct_300bp.bed \
../processed_data/basic_starrseq_results/minus_peaks_top10pct_300bp.bed \
../processed_data/basic_starrseq_results/plus_minus_peaks_top10pct_300bp.fasta \
../processed_data/basic_starrseq_results/negative_peaks_controls_top10pct_300bp.fasta
