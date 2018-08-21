import argparse
import string
from itertools import islice
import re

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


def tab_reader(filename):
	infile = open(filename)
	seqs = {}

	for line in infile:
		name, seq = line.strip().split('\t')
		seqs[name] = seq
		seqs[name] = seq

	return seqs


def add_stuffer(sequence, stuffer, tile_len):
	# store as tuple so we can easily insert the restriction site in between the
	# stuffer and tile
	stuffer = (stuffer[:(tile_len - len(sequence))], sequence)
	return stuffer


def parse_controls(filename):
	'''
	Extract promoters from CSV file. There are weird quotation marks and spaces
	so have to do some regex stuff
	'''

	infile = open(filename, 'r')

	syn_promoters = {}

	# read through header
	infile.readline()

	for line in infile.readlines():
		fields = line.strip().split(',')
		name = fields[0]
		seq = fields[9]

		# remove white space from seq
		seq = ''.join(seq.split())
		# remove quotations
		match = re.search('[ACGT]{1,}', seq)
		clean_seq = match.group(0)

		# remove quotations from name, always first and last characters
		clean_name = name.replace('\"', '')

		syn_promoters['pos_control_'+clean_name] = clean_seq

	return syn_promoters


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('predictions', help='uncertain predictions file, tab separated,\
		three columns: name, prediction, uncertainty. Assumes sorted in descending\
		order by uncertainty')
	parser.add_argument('lib', help='reference library, fasta format')
	parser.add_argument('lib_size', type=int, help="number of sequences")
	parser.add_argument('neg_controls', help='fasta file of negative controls (more than 200bp from TSS either strand)')
	parser.add_argument('pos_controls', help='fasta file of positive controls')
	parser.add_argument('output_name', help='name of output file')

	args = parser.parse_args()
	filename = args.predictions
	lib_size = args.lib_size
	neg_controls = tab_reader(args.neg_controls)
	pos_controls = parse_controls(args.pos_controls)
	output_name = args.output_name

	print "Reading in reference..."
	lib = fasta_reader(args.lib)
	print len(lib)

	print "Reading in predictions..."
	predictions = {}
	with open(filename) as infile:
		for i in range(lib_size+1):
			line = infile.readline()
			name, prediction, uncertainty = line.strip().split('\t')
			# remove '>'
			name = name[1:]
			predictions[name] = float(uncertainty)

	# skpp-160-F sk20mer-125313
	fwd_primer = 'TCATGCTGTGTCCATTAGCTCGAG'

	# skpp-285-R sk20mer-5049
	rev_primer = 'GCTAGCCAGAAGGTACGCTTTATG'

	tile_len = 150
	# T is easiest to synthesize, won't have secondary structure
	stuffer = 'T' * tile_len

	print "Writing library..."
	seqs = {x : lib[x] for x in predictions}

	# add controls to tiles
	seqs.update(neg_controls)
	
# positive controls are shorter, need to stuff
	for name, seq in pos_controls.items():
		if len(seq) < tile_len:
			seqs[name] = add_stuffer(seq, stuffer, tile_len)
		else:
			seqs[name] = seq

	xhoI = 'CTCGAG'
	with open(output_name, 'w') as outfile:
		for x in seqs:
			seq = seqs[x]
			if len(seq) == 2:
				stuffer, payload = seq
				# stuffed sequence, only use first 20bp of fwd primer, subtract 2bp from stuffer,
				# then add the full 6bp of RE (no overlap between primer and RE anymore)
				full_seq = fwd_primer[:20] + stuffer[:-2] + xhoI + payload + rev_primer
			else:
				full_seq = fwd_primer + seq + rev_primer
		
			reverse = best_A_content(full_seq)
			if reverse:
				x += '_flipped'
				full_seq = reverse_complement(full_seq)

			outfile.write(x + '\t' + full_seq + '\n')

	# with open(args.output_name, 'w') as outfile:
	# 	for x in predictions:
	# 		seq = lib.get(x, None)
	# 		if seq is not None:
	# 			full_seq = fwd_primer + seq + rev_primer
	# 			reverse = best_A_content(full_seq)
	# 			if reverse:
	# 				x += '_flipped'
	# 				full_seq = reverse_complement(full_seq)

	# 			outfile.write(x + '\t' + full_seq + '\n')
	# 		else:
	# 			missing += 1
	# 			continue

	# print missing, "missing sequences"



