import argparse
import string
from itertools import islice

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
	parser = argparse.ArgumentParser()
	parser.add_argument('predictions', help='uncertain predictions file, tab separated,\
		three columns: name, prediction, uncertainty. Assumes sorted in descending\
		order by uncertainty')
	parser.add_argument('lib', help='reference library, fasta format')
	parser.add_argument('lib_size', type=int, help="number of sequences")
	parser.add_argument('output_name', help='name of output file')

	args = parser.parse_args()
	filename = args.predictions
	lib = args.lib
	lib_size = args.lib_size

	print "Reading in reference..."
	lib = fasta_reader(lib)

	print "Reading in predictions..."
	predictions = {}
	with open(filename) as infile:
		for i in range(lib_size):
			line = infile.readline()
			name, prediction, uncertainty = line.strip().split('\t')
			# remove '>'
			name = name[1:]
			predictions[name] = float(uncertainty)

	# skpp-160-F sk20mer-125313
	fwd_primer = 'TCATGCTGTGTCCATTAGCTCGAG'

	# skpp-285-R sk20mer-5049
	rev_primer = 'GCTAGCCAGAAGGTACGCTTTATG'

	print "Writing library..."
	missing = 0
	with open(args.output_name, 'w') as outfile:
		for x in predictions:
			seq = lib.get(x, None)
			if seq is not None:
				full_seq = fwd_primer + seq + rev_primer
				reverse = best_A_content(full_seq)
				if reverse:
					x += '_flipped'
					full_seq = reverse_complement(full_seq)

				outfile.write(x + '\t' + full_seq + '\n')
			else:
				missing += 1
				continue

	print missing, "missing sequences"



