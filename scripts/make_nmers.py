# Make all kmers of length n and output in FASTA format, with header name also
# as the sequence itself

import itertools
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('n', help = 'length of kmer', type = int)
	parser.add_argument('output_name', help = 'name of output file')
	args = parser.parse_args()

	nts = ['A', 'C', 'G', 'T']

	with open(args.output_name, 'w') as outfile:
		for x in itertools.product(nts, repeat = args.n):
			nmer = ''.join(x)
			# fasta format, header then sequence
			outfile.write('>' + nmer + '\n')
			outfile.write(nmer + '\n')
