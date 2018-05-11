import argparse


def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
	
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	return rc


def library_reader(filename, primer_len, rev_complement=True, format = 'csv'):
	"""
	Read in .csv or tab file of library sequences. First column is sequence name,
	second column is sequence. Trim primer sequences.
	"""

	lib = {}

	with open(filename) as infile:
		# read through first line
		infile.readline()
		for line in infile:
			if format == 'csv':
				name, seq = line.strip().split(',')[:2]
			elif format == 'tab':
				name, seq = line.strip().split('\t')[:2]
			
			seq = seq.upper()
			if primer_len != 0:
				seq = seq[primer_len:-primer_len]

			if rev_complement:
				rc_seq = reverse_complement(seq)
				lib[rc_seq] = name

			lib[seq] = name

	return lib


def generate_all_mismatches(seq, name):
	nts = ['A', 'C', 'G', 'T']
	mismatches = []
	for i in range(len(seq)):
		for nt in nts:
			seq_list = list(seq)
			if nt != seq_list[i]:
				seq_list[i] = nt
				mismatch_name = name + '_pos' + str(i) + nt
				# mismatches[mismatch_name] = ''.join(seq_list)
				mismatches.append((mismatch_name, ''.join(seq_list)))
	return mismatches


def fasta_writer(seqs, outfile):
	for x in seqs:
		outfile.write(x[0]+'\n'+x[1]+'\n')


def tab_writer(seqs, outfile):
	for x in seqs:
		outfile.write(x[0]+'\t'+x[1]+'\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('lib_file', help='.csv or tab file of library')
	parser.add_argument('lib_type', help='please specify either csv or tab')
	parser.add_argument('primer_len', type = int, help='length of primers at \
		start and end of sequences in library file')
	parser.add_argument('output_name', help='name of output file')
	parser.add_argument('output_type', help='tab or fasta')

	args = parser.parse_args()

	lib_file = args.lib_file
	lib_type = args.lib_type
	primer_len = args.primer_len
	output_name = args.output_name
	output_type = args.output_type

	print "Reading in library reference..."
	lib = library_reader(
		filename=lib_file, 
		primer_len=primer_len,
		rev_complement=True, 
		format=lib_type)

	with open(output_name, 'w') as outfile:
		for x in lib:
			mismatches = generate_all_mismatches(x, lib[x])
			if output_type == 'tab':
				tab_writer(mismatches, outfile)
			elif output_type == 'fasta':
				fasta_writer(mismatches, outfile)
			else:
				raise Exception("Please specify tab or fasta for output type.")



