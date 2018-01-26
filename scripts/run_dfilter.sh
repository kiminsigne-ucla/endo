# ../bin/DFilter1.6/run_dfilter.sh -d=../rawdata/genome_frag/minus_RNA1.bed -c=../rawdata/genome_frag/minus_DNA1.bed -o=dfilter_peaks_minus_rep1.bed -type=bed -lpval=2 -ks=50 -bs=100 -dir -redund=100000
# ../bin/DFilter1.6/run_dfilter.sh -d=../rawdata/genome_frag/minus_RNA2.bed -c=../rawdata/genome_frag/minus_DNA2.bed -o=dfilter_peaks_minus_rep2.bed -type=bed -lpval=2 -ks=50 -bs=100 -dir -redund=100000

# ../bin/DFilter1.6/run_dfilter.sh -d=../rawdata/genome_frag/plus_RNA1.bed -c=../rawdata/genome_frag/plus_DNA1.bed -o=dfilter_peaks_plus_rep1.bed -type=bed -lpval=2 -ks=50 -bs=100 -dir -redund=100000
# ../bin/DFilter1.6/run_dfilter.sh -d=../rawdata/genome_frag/plus_RNA2.bed -c=../rawdata/genome_frag/plus_DNA2.bed -o=dfilter_peaks_plus_rep2.bed -type=bed -lpval=2 -ks=50 -bs=100 -dir -redund=100000

# combined replicates
../bin/DFilter1.6/run_dfilter.sh -d=../rawdata/genome_frag/Minus_RNA.bam -c=../rawdata/genome_frag/Minus_DNA.bam -o=dfilter_peaks_minus.bed -lpval=2 -ks=50 -bs=100 -dir -redund=100000 -f=bam -tagsize=auto -pm=10 -nonzero
../bin/DFilter1.6/run_dfilter.sh -d=../rawdata/genome_frag/Plus_RNA.bam -c=../rawdata/genome_frag/Plus_DNA.bam -o=dfilter_peaks_plus.bed -lpval=2 -ks=50 -bs=100 -dir -redund=100000 -f=bam -tagsize=auto -pm=10 -nonzero

mv *.bed ../processed_data/dfilter_results