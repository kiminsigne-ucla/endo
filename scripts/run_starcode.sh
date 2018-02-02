input=$1
topn=$2
output=$3

# sort by second column (weight from gkmSVM) in descending order
sort -k2 -n -r $input | 
# grab top k-mers
head -$topn |
# starcode requires raw one sequence per line
cut -f1 > tmp.txt

../bin/starcode/starcode -i tmp.txt -o $output --print-clusters -d 3
rm tmp.txt