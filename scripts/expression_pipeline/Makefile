MERGE_DEST=../../processed_data/endo_bc_map_merged_reads_pear
RAW_DATA=../../raw_data/endo_bc_map
DATA=../../processed_data/expression_pipeline
LIB=../../processed_data/lib_gen
DATE=$(shell date +%Y%m%d)

all: merge map


merge: $(MERGE_DEST)/endo_merged.fastq

map: $(DATA)/endo_mapping_barcode_statistics.txt

barcode_counts: $(DATA)/counts_rLP5_Endo2_DNA1.txt

active_tss_lib_gen: $(DATE)_active_tss_scrambled10_stride5.txt


$(MERGE_DEST)/pEndo2_84.paired.assembled.fastq: $(RAW_DATA)/pEndo2_84_R1.fastq $(RAW_DATA)/pEndo2_84_R2.fastq pearallel.sh
	./pearallel.sh $(RAW_DATA)/pEndo2_84 $(MERGE_DEST)


$(MERGE_DEST)/pEndo2_85.paired.assembled.fastq: $(RAW_DATA)/pEndo2_85_R1.fastq $(RAW_DATA)/pEndo2_85_R2.fastq pearallel.sh
	./pearallel.sh $(RAW_DATA)/pEndo2_85 $(MERGE_DEST)


$(MERGE_DEST)/endo_merged.fastq: $(MERGE_DEST)/pEndo2_84.paired.assembled.fastq $(MERGE_DEST)/pEndo2_85.paired.assembled.fastq
	@ cat $^ > $@


$(DATA)/endo_mapping_barcode_statistics.txt: $(MERGE_DEST)/endo_merged.fastq
	@ python barcode_mapping.py $< fastq ../../ref/endo_lib_2016_controls_clean.txt \
	tab 150 24 end 20 $(DATA)/endo_mapping --controls ../../ref/synthetic_promoter_pos_controls_formatted.txt \
	--controls_primer_len 0


# use one FASTQ file to trigger barcode counting for all FASTQ files
$(DATA)/counts_rLP5_Endo2_DNA1.txt: ../../raw_data/endo_DNA_RNA/rLP5_Endo2_DNA1.fastq
	@ ./barcode_counting.sh $(RAW_DATA)
	@ mv counts* $(DATA)

# calculate expression
$(DATA)/rLP5_Endo2_expression.txt: $(DATA)/counts_rLP5_Endo2_DNA1.txt \
	$(DATA)/endo_mapping_barcode_statistics.txt
	@ Rscript bc_counts2expression.R

# threshold active TSS
$(DATA)/tss_positives.bed: $(DATA)/rLP5_Endo2_expression.txt
	@ Rscript threshold_tss.R

# convert bed to fasta	
$(DATA)/tss_positives.fasta: $(DATA)/tss_positives.bed
	@ bedtools getfasta -fi ../../ref/U00096.2.fasta -bed $(DATA)/tss_positives.bed \
	-fo $(DATA)/tss_positives.fasta -name -s

# generate active TSS scramble library
$(LIB)/$(DATE)_active_tss_scrambled10_stride5.txt: $(DATA)/tss_positives.fasta
	@ echo "Generating scramble library..."
	@ python scramble_tiling.py $< $(LIB)/neg_controls_no_primers.txt \
	../../ref/synthetic_promoter_pos_controls.csv 10 5 $@
