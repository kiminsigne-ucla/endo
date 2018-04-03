DATE=`date +%Y%m%d`
python peak_tiling.py ../../processed_data/custom_peaks/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.fasta \
neg_controls_no_primers.txt synthetic_promoter_pos_controls.csv \
10 150 ${DATE}_lb_peak_tile_lib.txt