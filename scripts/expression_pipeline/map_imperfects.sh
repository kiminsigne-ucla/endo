DATA=../../processed_data/expression_pipeline/imperfects

# echo "Generating imperfects..."
# python generate_imperfects.py ../../ref/endo_lib_2016_controls_clean.txt tab \
# 24 ${DATA}/endo_lib_imperfects.fasta fasta

# # set cutoff, otherwise randomly choosing sequences in bootstrapping from large
# # library takes way too long
# python barcode_mapping.py ../../processed_data/endo_bc_map_merged_reads_pear/endo_merged.fastq \
# fastq endo_lib_imperfects.txt tab 150 0 end 20 ${DATA}/endo_mapping_imperfect --cutoff 70

# Rscript bc_counts2expression.R ../../processed_data/expression_pipeline \
#  ${DATA}/endo_mapping_imperfect_barcode_statistics.txt \
#  ${DATA}/endo_mapping_imperfect_variant_statistics.txt \
#  ../../ref/endo_lib_2016_controls_clean.txt ${DATA}/endo_imperfect_expression.txt

# predict with uncertainty
echo "Predicting..."
touch endo_imperfect_uncertain_predictions.txt
for i in `seq 1 100000 200001`
do
	stop=$((i+ 100001))
	echo "Predicting $i to $stop"
	sed -n $i,${stop}p ${DATA}/endo_lib_imperfects.fasta > chunk.fasta
	python ../nn_regression_with_uncertainty.py \
	../hyperparameter_best_uncertainty_trained_model.arch.json \
	../hyperparameter_best_uncertainty_trained_model.weights.h5 \
	chunk.fasta 100 tmp.txt
	cat tmp.txt >> endo_imperfect_uncertain_predictions.txt
done

# 3256500
