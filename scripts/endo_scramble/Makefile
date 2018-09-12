MERGE_DEST=../../processed_data/endo_scramble
RAW_DATA=../../raw_data/endo_scramble
DATA=../../processed_data/endo_scramble


merge: $(MERGE_DEST)/endo_scramble_merged.fastq

map: $(DATA)/endo_scramble_bc_map_consensus.txt


################################ Read merging ################################
$(MERGE_DEST)/pEndoScramble_80.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndoScramble_80_S1_R1_001.fastq \
	$(RAW_DATA)/bc_map/pEndoScramble_80_S1_R2_001.fastq ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndoScramble_80 $(MERGE_DEST)


$(MERGE_DEST)/pEndoScramble_81.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndoScramble_81_S2_R1_001.fastq \
	$(RAW_DATA)/bc_map/pEndoScramble_81_S2_R2_001.fastq ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndoScramble_81 $(MERGE_DEST)


$(MERGE_DEST)/endo_scramble_merged.fastq: $(MERGE_DEST)/pEndoScramble_80.paired.assembled.fastq $(MERGE_DEST)/pEndoScramble_81.paired.assembled.fastq
	@ cat $^ > $@


################################ Mapping ####################################
$(DATA)/endo_scramble_bc_map_consensus.txt: $(MERGE_DEST)/endo_scramble_merged.fastq
	@ python ../bc_map_consensus.py $< ../../ref/20180507_active_tss_scrambled10_stride5_trimmed.fasta \
	--bc-start -20 --bc-length 20 --proc 30 --verbose --bbmap-procs 30 \
	--var_start 1 --var_length 150 --min-reads 3 --bad-bcs $(DATA)/bc_map_consensus_bad_bcs.txt > $@
