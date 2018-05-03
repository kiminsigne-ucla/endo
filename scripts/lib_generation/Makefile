PEAK_DEST=../../processed_data/custom_peaks

pileup: plus_frag_pileup.wig

call_peaks: $(PEAK_DEST)/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.fasta


# fragment expression pileup - get expression of each fragment at every position, this will create both strands
plus_frag_pileup.wig: ../../processed_data/U00096.2_frag-rLP5_LB_expression.txt
	@ python ../frag_expression_pileup.py $< frag_pileup.txt 50


# call peaks based on fragment expression pileup
$(PEAK_DEST)/U00096.2_plus_called_peaks_threshold1.1_merge20_min60.bed: plus_frag_pileup.wig
	python ../call_peaks.py $< 1.1 20 60 + $@


$(PEAK_DEST)/U00096.2_minus_called_peaks_threshold1.1_merge20_min60.bed: minus_frag_pileup.wig
	python ../call_peaks.py $< 1.1 20 60 - $@


$(PEAK_DEST)/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.bed: \
	$(PEAK_DEST)/U00096.2_plus_called_peaks_threshold1.1_merge20_min60.bed \
	$(PEAK_DEST)/U00096.2_minus_called_peaks_threshold1.1_merge20_min60.bed
	 @ echo "Combine peaks..."
	 @ cat $^ > $@


$(PEAK_DEST)/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.fasta: \
	$(PEAK_DEST)/U00096.2_plus_minus_called_peaks_threshold1.1_merge20_min60.bed
	@ echo "Converting to FASTA..."
	@ bedtools getfasta -fi ../ref/U00096.2.fasta -bed  $< -fo $@ -name -s





