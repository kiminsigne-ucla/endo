num_layer=$1
python dragonn_hyperparameter_tuning.py ../processed_data/tss_positives.fasta ../processed_data/tss_negatives.fasta 150 $num_layer 0.2 0.2 50 > hyperparameter_log_layer$num_layer.txt

