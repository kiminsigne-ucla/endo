num_layer=$1
min_filter=$2
max_filter=$3
DATE=`date +%Y%m%d`
python -u dragonn_hyperparameter_tuning.py ../../processed_data/endo_tss/lb/tss_positives.fasta \
../../processed_data/endo_tss/lb/tss_negatives.fasta 150 $num_layer $min_filter $max_filter \
0.2 0.2 50 > ${DATE}_endo_tss_classifier.log

