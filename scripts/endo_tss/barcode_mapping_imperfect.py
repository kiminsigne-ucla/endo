"""
This script selects barcodes that will be used for downstream analysis
in RNA-seq. The perfect length requirement is not used and the Levenshtein
distance is set to the 1st percentile, not (arbitrarily) at 5
"""

import sys
import os
import itertools
from collections import Counter, defaultdict
import re
import Levenshtein
import numpy
import random
import argparse
import subprocess


def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
	
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	return rc


def fastq_reader(filename):
	"""
	Extract sequences from FASTQ file. Each sequence takes up four lines,
	with the sequence being the second line.
	"""

	with open(filename) as infile:
		# grab lines starting at second line (1, 0-indexing), go until end of
		# file (stop = None), skip 4 lines at a time 
		line_slice = itertools.islice(infile, 1, None, 4)
		
		# convert iterator to list, strip new lines
		sequences = []
		for line in line_slice:
			sequences.append(line.strip())

	return sequences


def fasta_reader(filename):

	seqs = {}

	with open(filename) as infile:
		seq = ''
		for line in infile:
			if line.startswith('>'):
				if len(seq) != 0:
					seqs[name] = seq
				
				name = line.strip()[1:] # remove leading '>'
				seq = ''
			else:
				seq += line.strip()

		# catch last sequence
		if len(seq) != 0:
			seqs[name] = seq

	return seqs


def get_wc(filename):
	p = subprocess.Popen(['wc', '-l', filename], stdout = subprocess.PIPE,
												 stderr = subprocess.PIPE)
	out, err = p.communicate()
	lines = float(out.split()[0])
	return lines


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


def find_imperfect_reads(reads, lib, var_len, bc_loc):
	"""
	Only take reads that have one substitution from reference
	"""
	# return imperfect_reads
	imperfect_reads = []
	for i in range(len(reads)):
		read = reads[i]
		if i % 1000000 == 0:
			print i
		if bc_loc == 'start':
			variant = read[-var_len:]
		if bc_loc == 'end':
			variant = read[:var_len]

		# calculate distance between variant and each library member
		for x in lib:
			if Levenshtein.distance(variant, x) == 1:
				imperfect_reads.append(read)

	return imperfect_reads


def mapping(barcodes, perfect_reads, reads, bc_loc, bc_len, var_len):

	variant_map = defaultdict(list)
	barcode_map = defaultdict(list)

	# for each barcode that passes filters, look at all reads and see what
	# it maps to
	for read in reads:
		if bc_loc == 'start':
			barcode = read[:bc_len]
		elif bc_loc == 'end':
			barcode = read[-bc_len:]

		# if barcode in filtered set
		if barcodes.get(barcode, 0) > 0:
			barcode_map[barcode].append(read)

	# in the perfect reads, keep track of how many barcodes go to each variant
	for read in perfect_reads:
		if bc_loc == 'start':
			barcode = read[:bc_len]
			variant = read[-var_len:]
		elif bc_loc == 'end':
			variant = read[:var_len]
			barcode = read[-bc_len:]

		variant_map[variant].append(barcode)

	return [variant_map, barcode_map]


def bootstrap_levenshtein(lib, n):
	"""
	This function calculates a reference Levenshtein distribution. It randomly
	picks two sequences from the reference sequences and calculates the distance
	to get a measure of how similar the library is.
	"""

	distances = []
	# bootstrap n times
	for i in range(0, n):
		# randomly grab two sequences with replacement
		string1 = random.choice(lib.keys())
		string2 = random.choice(lib.keys())

		distances.append(Levenshtein.distance(string1, string2))
	
	# take cutoff at 1% percentile
	cutoff = numpy.percentile(distances, 1)
	 
	return cutoff


def filter_barcodes(barcode_map, cutoff, bc_loc, bc_len, var_len, name, lib):
	'''
	For each barcode, calculate the Levenshtein distance between its reads
	and if it is below cutoff (aka barcode maps to similar reads, no cross-talk)
	then keep this barcode
	'''

	final_barcodes = []
	covered_sequences = set()
	lib = set(lib.keys()) # for faster lookup
	if controls:
		lib = lib.update(set(controls.keys()))

	all_dist = []

	outfile = open(name, 'w')
	headers = ['barcode', 'num_unique', 'num_reads', 'num_reads_most_common', 
	'most_common']
	outfile.write('\t'.join(headers)+'\n')

	for barcode in barcode_map:
		reads = barcode_map[barcode]
		# # trim off barcode (20), RE site (8), primers (15) and last 15
		# trimmed = [read[43:-15] for read in reads]
		if bc_loc == 'start':
			trimmed = [read[-var_len:] for read in reads]
		if bc_loc == 'end':
			trimmed = [read[:var_len] for read in reads]
		# grab most common read as reference
		most_common = Counter(trimmed).most_common(1)[0][0]

		# calculate distance between each read and reference
		distances = [Levenshtein.distance(most_common, read) for read in set(trimmed)]
		all_dist.append(max(distances))
		# max distance for set of reads belonging to a barcode must be below cutoff
		if max(distances) < cutoff:
			num_unique = len(set(trimmed))
			num_reads = len(trimmed)
			num_reads_most_common = Counter(trimmed).most_common(1)[0][1]
			if bc_loc == 'start':
				most_common = Counter([read[-var_len:] for read in reads]).most_common(1)[0][0]
			if bc_loc == 'end':
				most_common = Counter([read[:var_len] for read in reads]).most_common(1)[0][0]
			if most_common in lib:
				# only accept barcode if the most common read is in the library
				final_barcodes.append(barcode)
				# is_reference = 1
				covered_sequences.add(most_common)
				info = [barcode, num_unique, num_reads, num_reads_most_common, most_common]
				info = map(str, info)
				outfile.write('\t'.join(info)+'\n')

	outfile.close()

	coverage = len(covered_sequences)/ float((len(lib)/2.0))
	print "Percent of library represented by final barcodes:", coverage

	return final_barcodes


def write_variant_results(variant_map, name, final_barcodes, lib):
	outfile = open(name, 'w')
	fields = ['variant', 'name', 'num_barcodes', 'num_barcodes_unique', 'barcodes']
	outfile.write('\t'.join(fields)+'\n')

	final_barcodes = set(final_barcodes)

	for variant in variant_map:

		barcodes = variant_map[variant]
		# only keep those in final barcodes
		barcodes = [barcode for barcode in barcodes if barcode in final_barcodes]
		num_barcodes = len(barcodes)
		uniq_bcs = set(barcodes)
		num_unique = len(uniq_bcs)
		ref_name = lib.get(variant, None)
		if not ref_name:
			continue
		# ref_name = lib[variant]
		info = [variant, ref_name, num_barcodes, num_unique]
		info = map(str, info)
		info.append(','.join(uniq_bcs))
		outfile.write('\t'.join(info)+'\n')

	outfile.close()


def check_args(args):
	if args.lib_type != 'csv' and args.lib_type != 'tab':
		raise ValueError('Please provide format type of library, either csv or tab') 

	if args.bc_loc != 'start' and args.bc_loc != 'end':
		raise ValueError('Please specify location of barcode, either start or end')


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Map barcodes to sequence')
	parser.add_argument('reads_file', help='sequence reads, accepted formats \
		FASTQ, FASTA or plain-text, one read per line')
	parser.add_argument('file_type', help='sequence read format, \
		specify fastq/fasta/txt')
	parser.add_argument('lib_file', help='.csv or tab file of library')
	parser.add_argument('lib_type', help='please specify either csv or tab')
	parser.add_argument('var_len',type=int, help='length of variants')
	parser.add_argument('primer_len', type = int, help='length of primers at \
		start and end of sequences in library file')
	parser.add_argument('bc_loc', help='barcode location, specify start or end')
	parser.add_argument('bc_len', type=int, help='length of barcode')
	parser.add_argument('output_prefix', help='Name of output file prefix')
	parser.add_argument('--cutoff', type=int, help='user defined Levenshtein cutoff, \
		if not given then empirically determined (bootstrapped')
	args = parser.parse_args()

	check_args(args)

	reads_file = args.reads_file
	file_type = args.file_type
	lib_file = args.lib_file
	lib_type = args.lib_type
	var_len = args.var_len
	primer_len = args.primer_len
	bc_loc = args.bc_loc
	bc_len = args.bc_len
	output_prefix = args.output_prefix

	if file_type == 'fastq':
		reads = fastq_reader(reads_file)
	elif file_type == 'fasta':
		# returns dictionary, just need reads (values) and not names (keys)
		reads = fasta_reader(reads_file).values()
	elif file_type == 'txt':
		reads = [line.strip() for line in open(reads_file)]
	else:
		raise Exception("Please specify either fastq, fasta or txt (raw one read per line)")

	print "Number of reads:", len(reads)

	print "Reading in library reference..."
	lib = library_reader(
		filename=lib_file, 
		primer_len=primer_len,
		rev_complement=True, 
		format=lib_type)

	print "Extracting perfect reads..."
	perfect_reads = find_perfect_reads(
		reads=reads, 
		lib=lib, 
		var_len=var_len, 
		bc_loc=bc_loc,
		controls=controls)	

	print "Percent perfect:", len(perfect_reads) / float(len(reads))

	# if controls:
	# 	perfect_control_reads = find_perfect_controls(reads, controls)
	# 	print "Percent perfect controls:", len(perfect_control_reads) / float(len(reads))
	# 	perfect_reads.extend(perfect_control_reads)

	# grab barcodes that map to a perfect sequence
	if bc_loc == 'start':
		barcodes = [read[:bc_len] for read in perfect_reads]
	elif bc_loc == 'end':
		barcodes = [read[-bc_len:] for read in perfect_reads]

	print "Number of unique barcodes for perfect reads: ", len(set(barcodes))
	print "Filter by barcode frequency..."
	
	# Count barcodes 
	barcode_counts = dict(Counter(barcodes))

	# Throw out barcodes that appear 1 or 2 times, sequencing errors
	barcodes_clean = {x : barcode_counts[x] for x in barcode_counts if barcode_counts[x] > 2}
	print "Number of barcodes > 2:", len(barcodes_clean)

	# barcode_cutoff = args.bc_cutoff
	# above_cutoff = {x : barcodes_clean[x] for x in barcodes_clean if barcodes_clean[x] >= barcode_cutoff}
	
	# print "Number of barcodes above cutoff:", len(above_cutoff)

	print "Mapping..."

	variant_map, barcode_map = mapping(
		barcodes_clean, 
		perfect_reads, 
		reads, 
		bc_loc, 
		bc_len, 
		var_len)

	# bootstrap reference sequences to get a reference Levenshtein distribution 
	# to determine cutoff
	if args.cutoff:
		cutoff = args.cutoff
	else:
		print "Bootstrapping reference sequences to obtain cutoff...", 
		cutoff = bootstrap_levenshtein(lib, 10000)
		print "cutoff is Levenshtein distance ", cutoff

	print "Filtering and writing results..."
	final_barcodes = filter_barcodes(
		barcode_map, 
		cutoff,
		bc_loc,
		bc_len,
		var_len,
		output_prefix + '_barcode_statistics.txt', 
		lib)

	print "Number of final barcodes: ", len(final_barcodes)
	# write final barcodes to file
	outfile = open(output_prefix + '_mapped_barcodes.txt', 'w')
	for barcode in final_barcodes:
		outfile.write(barcode+'\n')
	outfile.close()

	# write variant results
	write_variant_results(
		variant_map, 
		output_prefix + '_variant_statistics.txt', 
		final_barcodes, 
		lib, 
		controls)

