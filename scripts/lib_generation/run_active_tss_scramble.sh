scramble_len=$1
stride_len=$2
DATE=`date +%Y%m%d`
python scramble_tiling.py ../../processed_data/tss_positives.fasta $scramble_len $stride_len \
${DATE}_active_tss_scrambled${scramble_len}_stride${stride_len}.txt