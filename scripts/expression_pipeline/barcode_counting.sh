echo "Extracting & counting barcodes..."

for i in rLP5*.fastq; do awk 'NR%4==2' $i |\
	 cut -c 1-20 | rev | tr ACGT TGCA |\
	 sort -T ./ --parallel=20 | uniq -c > counts_${i/.fastq/}.txt; done
