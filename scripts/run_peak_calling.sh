# forward strand
echo "Calling peaks for forward strand..."
python call_peaks.py ../processed_data/expn_wig/U00096.2_plus_expression.wig 1.1 20 60 + ../processed_data/custom_peaks/plus_called_peaks_threshold1.1_merge20_min60.bed
# reverse strand
echo "Calling peaks for reverse strand..."
python call_peaks.py ../processed_data/expn_wig/U00096.2_minus_expression.wig 1.1 20 60 - ../processed_data/custom_peaks/minus_called_peaks_threshold1.1_merge20_min60.bed
# combine
echo "Combine peaks..."
cat ../processed_data/custom_peaks/plus_called_peaks_threshold1.1_merge20_min60.bed \
../processed_data/custom_peaks/minus_called_peaks_threshold1.1_merge20_min60.bed > \
../processed_data/custom_peaks/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.bed

# convert to FASTA
echo "Converting to FASTA..."
bedtools getfasta -fi ../ref/U00096.2.fasta -bed ../processed_data/custom_peaks/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.bed \
-fo ../processed_data/custom_peaks/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.fasta -name -s