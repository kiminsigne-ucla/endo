MERGE_DEST=../../processed_data/merged_reads_pear


all: merge map


merge: $(MERGE_DEST)/endo_merged.fastq

map: ../../processed_data/endo_bc_mapping.txt	


$(MERGE_DEST)/pEndo2_84.paired.assembled.fastq: pEndo2_84_R1.fastq pEndo2_84_R2.fastq pearallel.sh
	./pearallel.sh pEndo2_84 $(MERGE_DEST)


$(MERGE_DEST)/pEndo2_85.paired.assembled.fastq: pEndo2_85_R1.fastq pEndo2_85_R2.fastq pearallel.sh
	./pearallel.sh pEndo2_85 $(MERGE_DEST)


$(MERGE_DEST)/endo_merged.fastq: $(MERGE_DEST)/pEndo2_84.paired.assembled.fastq $(MERGE_DEST)/pEndo2_85.paired.assembled.fastq
	@ cat $^ > $@


../../processed_data/endo_bc_mapping.txt: $(MERGE_DEST)/endo_merged.fastq
	@ python barcode_mapping.py $< fastq ../../ref/endo_lib_2016_controls_clean.txt \
	tab 150 24 end 20 $@ --cutoff 5