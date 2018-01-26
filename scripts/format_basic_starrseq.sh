plus_bed=$1
minus_bed=$2
output=$3
negative=$4

# combine strands
cat $plus_bed $minus_bed > plus_minus.bed

FASTA="../ref/Escherichia_coli_K-12_MG1655.fasta"
FAI="../ref/Escherichia_coli_K-12_MG1655.fasta.fai"

# Get fasta sequences of peak calls from .bed file
bedtools getfasta -fi $FASTA -bed plus_minus.bed -fo $output -s

# Generate bed file containing sequences that are of the same size as the peak calls, 
# but with no overlap  (Negative Controls), then convert to fasta
bedtools shuffle -i plus_minus.bed -g $FAI -seed 123123123 -noOverlapping -excl plus_minus.bed | 
bedtools getfasta -fi $FASTA -bed - -fo $negative -s