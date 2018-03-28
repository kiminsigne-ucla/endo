import argparse
from itertools import islice
import string
import numpy as np
import random


def fasta_reader(filename):
	"""
	Input: str, name of file
	Output: dictionary, key = header, value = sequence
	"""

	seqs = {}
	with open(filename) as infile:
		entry = list(islice(infile, 2))
		while entry:
			# grab two lines from file at a time, strip \n
			header, seq = map(string.strip, entry)
			# strip '>' from header, first character
			seqs[header[1:]] = seq.upper()
			entry = list(islice(infile, 2))

		return seqs


def tab_reader(filename):
	infile = open(filename)
	seqs = {}

	for line in infile:
		name, seq = line.strip().split('\t')
		seqs[name] = seq
		seqs[name] = seq

	return seqs


def best_A_content(oligo):
	'''
	Choose the strand with the lowest A content because A's are harder to
	synthesize. Return true if sequence needs to be reverse complemented
	'''
	rc_oligo = reverse_complement(oligo)

	oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
	rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

	if oligo_As < rc_As:
		return False
	else:
		return True


def reverse_complement(sequence):
    """Return the reverse complement of the supplied sequence string """ 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', '.':'.', '*':'*'} 
    #Reverse, convert to uppercase
    sequence = sequence[::-1].upper()
    #Complement
    letters = list(sequence) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters) 


if __name__ == '__main__':
	parser = argparse.ArgumentParser("script to generate scrambled regions tiling input sequences")
	parser.add_argument('sequences', help='filename of input sequences, FASTA format')
	# parser.add_argument('neg_sequences', help='file of negative controls (more than 200bp from TSS either strand), used to pad shorter sequences, tab format')
	parser.add_argument('scramble_len', type=int, help='length of scrambled segments')
	parser.add_argument('output_name', help='name of output file')

	args = parser.parse_args()
	sequences = fasta_reader(args.sequences)
	# stuffer_sequences = tab_reader(args.neg_sequences)
	scramble_len = args.scramble_len
	output_name = args.output_name

	# skpp-100-F sk20mer-20955
	fwd_primer = 'ACCTGTAATTCCAAGCGTCTCGAG'
	# skpp-155-R sk20mer-121927
	rev_primer = 'GCTAGCGGTGTTTAGTTAGCATCC'

	tiles = {}

	for name, seq in sequences.items():
		for i in range(0, len(seq), scramble_len)
			scrambled = list(seq[i:i+scramble_len])
			random.shuffle(scrambled)
			tile = seq[:i] + ''.join(scrambled).lower() + seq[i+scramble_len:]
			# print ''.join(scrambled).lower()
			tile_name = name + '_scrambled' + str(i) + '-' + str(i+scramble_len)
			tiles[tile_name] = tile

	with open(output_name, 'w') as outfile:
		for tile_name in sorted(tiles.keys()):
			tile = tiles[tile_name]
			full_tile = fwd_primer + tile + rev_primer
			reverse = best_A_content(full_tile)
			if reverse:
				tile_name += '_flipped'
				full_tile = reverse_complement(full_tile)
			outfile.write(tile_name + '\t' + tile + '\n')


