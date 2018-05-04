MERGE_DEST=../../processed_data/merged_reads_pear


all: merge map


merge: $(MERGE_DEST)/endo_merged.fastq

map: endo_mapping_barcode_statistics.txt

barcode_counts: final_bc/counts_rLP5_Endo2_DNA1.txt


$(MERGE_DEST)/pEndo2_84.paired.assembled.fastq: pEndo2_84_R1.fastq pEndo2_84_R2.fastq pearallel.sh
	./pearallel.sh pEndo2_84 $(MERGE_DEST)


$(MERGE_DEST)/pEndo2_85.paired.assembled.fastq: pEndo2_85_R1.fastq pEndo2_85_R2.fastq pearallel.sh
	./pearallel.sh pEndo2_85 $(MERGE_DEST)


$(MERGE_DEST)/endo_merged.fastq: $(MERGE_DEST)/pEndo2_84.paired.assembled.fastq $(MERGE_DEST)/pEndo2_85.paired.assembled.fastq
	@ cat $^ > $@


endo_mapping_barcode_statistics.txt: $(MERGE_DEST)/endo_merged.fastq
	@ python barcode_mapping.py $< fastq ../../ref/endo_lib_2016_controls_clean.txt \
	tab 150 24 end 20 endo_mapping --controls ../../ref/synthetic_promoter_pos_controls_formatted.txt \
	--controls_primer_len 0


# use one FASTQ file to trigger barcode counting for all FASTQ files
counts_rLP5_Endo2_DNA1.txt: rLP5_Endo2_DNA1.fastq
	@ ./barcode_counting.sh