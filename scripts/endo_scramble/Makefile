MERGE_DEST=../../processed_data/endo_scramble
RAW_DATA=../../raw_data/endo_scramble
DATA=../../processed_data/endo_scramble


merge: $(MERGE_DEST)/endo_scramble_merged.fastq

map: $(DATA)/endo_scramble_mapping_barcode_statistics.txt


$(MERGE_DEST)/pEndoScramble_80.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndoScramble_80_S1_R1_001.fastq \
	$(RAW_DATA)/bc_map/pEndoScramble_80_S1_R2_001.fastq ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndoScramble_80 $(MERGE_DEST)


$(MERGE_DEST)/pEndoScramble_81.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndoScramble_81_S2_R1_001.fastq \
	$(RAW_DATA)/bc_map/pEndoScramble_81_S2_R2_001.fastq ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndoScramble_81 $(MERGE_DEST)


$(MERGE_DEST)endo_scramble_merged.fastq: $(MERGE_DEST)/pEndoScramble_80.paired.assembled.fastq $(MERGE_DEST)/pEndoScramble_81.paired.assembled.fastq
	@ cat $^ > $@


$(DATA)/endo_scramble_mapping_barcode_statistics.txt: $(MERGE_DEST)/endo_scramble_merged.fastq
	@ python -u ../barcode_mapping.py $< fastq ../../ref/20180507_active_tss_scrambled10_stride5.txt \
	tab 150 24 end 20 $(DATA)/endo_scramble_mapping --controls ../../ref/synthetic_promoter_pos_controls_formatted.txt \
	--controls_primer_len 0 > $(DATA)/mapping_log_endo_scramble.txt
