MERGE_DEST=../../processed_data/peak_tile
RAW_DATA=../../raw_data/peak_tile
DATA=../../processed_data/peak_tile


merge: $(MERGE_DEST)/peak_tile_merged.fastq

map: $(DATA)/peak_tile_mapping_barcode_statistics.txt


$(MERGE_DEST)/pPeakTiling_82.paired.assembled.fastq: $(RAW_DATA)/bc_map/pPeakTiling_82_S3_R1_001.fastq \
	$(RAW_DATA)/bc_map/pPeakTiling_82_S3_R2_001.fastq ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pPeakTiling_82 $(MERGE_DEST)


$(MERGE_DEST)/pPeakTiling_83.paired.assembled.fastq: $(RAW_DATA)/bc_map/pPeakTiling_83_S4_R1_001.fastq \
	$(RAW_DATA)/bc_map/pPeakTiling_83_S4_R2_001.fastq ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pPeakTiling_83 $(MERGE_DEST)


$(MERGE_DEST)/peak_tile_merged.fastq: $(MERGE_DEST)/pPeakTiling_82.paired.assembled.fastq $(MERGE_DEST)/pPeakTiling_83.paired.assembled.fastq
	@ cat $^ > $@


$(DATA)/peak_tile_mapping_barcode_statistics.txt: $(MERGE_DEST)/peak_tile_merged.fastq
	@ python ../barcode_mapping.py $< fastq ../../ref/20180508_lb_peak_tile_lib.txt \
	tab 150 24 end 20 $(DATA)/peak_tile_mapping --controls ../../ref/synthetic_promoter_pos_controls_formatted.txt \
	--controls_primer_len 0 > $(DATA)/mapping_log_peak_tile.txt
