plus_bed=$1
minus_bed=$2
output=$3
negative=$4

# # convert to appropriate columns for get fasta
# # chr, start, end, name, 1, strand
# # Convert all plus peak calls to .bed format
# awk '{print $1, $2, $3, $2$3"+", "1", "+"}' $plus_bed | uniq | sed 's/  */\t/g' > a.txt
# # Same for minus peak calls
# awk '{print $1, $2, $3, $2$3"-", "1", "-"}' $minus_bed | uniq | sed 's/  */\t/g'  > b.txt 
# # Combine + and - strand
# cat a.txt b.txt > stranded_150bp_.95qt_BasicStarrPeaks.bed
# rm -f a.txt b.txt

# combine strands
cat $plus_bed $minus_bed > plus_minus.bed

FASTA="../../ref/Escherichia_coli_K-12_MG1655.fasta"
FAI="../../ref/Escherichia_coli_K-12_MG1655.fasta.fai"
# BED="stranded_150bp_.95qt_BasicStarrPeaks.bed"

# Get fasta sequences of peak calls from .bed file
bedtools getfasta -fi $FASTA -bed plus_minus.bed -fo $output -s

# Generate bed file containing sequences that are of the same size as the peak calls, but with no overlap  (Negative Controls)
bedtools shuffle -i plus_minus.bed -g $FAI-seed 123123123 -noOverlapping -excl plus_minus.bed > $negative.bed

# Extract ‘Negative Sequences’ from Fasta
bedtools getfasta -fi $FASTA -bed $negative.bed -fo $negative.fasta -s