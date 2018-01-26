# minus strand peak calling
macs2 callpeak -t ../rawdata/genome_frag/Minus_RNA.bam -c ../rawdata/genome_frag/Minus_DNA.bam -f BAM -g 9.2e6 \
--keep-dup all -n minus_frag --nomodel --broad


# plus strand peak calling
macs2 callpeak -t ../rawdata/genome_frag/Plus_RNA.bam -c ../rawdata/genome_frag/Plus_DNA.bam -f BAM -g 9.2e6 \
--keep-dup all -n plus_frag --nomodel --broad