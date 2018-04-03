scramble_len=$1
stride_len=$2
DATE=`date +%Y%m%d`
python scramble_tiling.py ../../processed_data/tss_positives.fasta neg_controls_no_primers.txt \
synthetic_promoter_pos_controls.csv $scramble_len $stride_len \
${DATE}_active_tss_scrambled${scramble_len}_stride${stride_len}.txt