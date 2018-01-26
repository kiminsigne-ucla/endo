# run peak calling in R
Rscript run_basic_starrseq.R

# format results
./format_basic_starrseq.sh ../processed_data/basic_starrseq_results/plus_peaks_100bp.bed \
../processed_data/basic_starrseq_results/minus_peaks_100bp.bed \
../processed_data/basic_starrseq_results/plus_minus_peaks_100bp.fasta \
../processed_data/basic_starrseq_results/negative_peak_controls.fasta
