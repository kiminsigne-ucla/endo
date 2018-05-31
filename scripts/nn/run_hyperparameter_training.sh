num_layer=$1
min_filter=$2
max_filter=$3
DATE=`date +%Y%m%d`
python -u dragonn_hyperparameter_tuning.py ../../processed_data/expression_pipeline/tss_positives.fasta \
../../processed_data/expression_pipeline/tss_negatives.fasta 150 $num_layer $min_filter $max_filter \
0.2 0.2 50 > ../logs/${DATE}_hyperparameter_log_layer${num_layer}_minfilter${min_filter}_maxfilter${max_filter}.txt

