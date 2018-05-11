DATA=../../processed_data/expression_pipeline/imperfects

# python generate_imperfects.py ../../ref/endo_lib_2016_controls_clean.txt tab \
# 24 ${DATA}/endo_lib_imperfects.txt

# # set cutoff, otherwise randomly choosing sequences in bootstrapping from large
# # library takes way too long
# python barcode_mapping.py ../../processed_data/endo_bc_map_merged_reads_pear/endo_merged.fastq \
# fastq endo_lib_imperfects.txt tab 150 0 end 20 ${DATA}/endo_mapping_imperfect --cutoff 70

Rscript bc_counts2expression.R ../../processed_data/expression_pipeline \
 ${DATA}/endo_mapping_imperfect_barcode_statistics.txt \
 ${DATA}/endo_mapping_imperfect_variant_statistics.txt \
 ../../ref/endo_lib_2016_controls_clean.txt ${DATA}/endo_imperfect_expression.txt
