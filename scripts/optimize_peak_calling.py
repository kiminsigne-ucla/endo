import argparse
import os
import call_peaks
import numpy as np


def count_overlap(file1, file2):
	# by default, if an overlap is found, reports shared interval between two
	# overlapping features
	cmd = 'bedtools intersect -a ' + file1 + ' -b ' + file2 + ' > overlap.bed'
	os.system(cmd)
	# count overlap
	num_overlap = sum(1 for line in open('overlap.bed'))
	return num_overlap


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('pos_file', help='bed file of positive TSSs')
	parser.add_argument('neg_file', help='bed file of negative TSSs')
	parser.add_argument('plus_wig', help='wig file for plus strand')
	parser.add_argument('minus_wig', help='wig file for minus strand')
	parser.add_argument('min_threshold', type=float, help='minimum threshold')
	parser.add_argument('max_threshold', type=float, help='maximum threshold')
	parser.add_argument('threshold_step', type=float, help='step size for threshold')

	args = parser.parse_args()
	pos_file = args.pos_file
	neg_file = args.neg_file
	min_threshold = args.min_threshold
	max_threshold = args.max_threshold
	step = args.threshold_step

	# separate positive and negative files by strand
	cmd = 'awk \'{if ($6 == \"+\") print $0}\' ' + pos_file + ' > pos_file_plus.bed'
	os.system(cmd)

	cmd = 'awk \'{if ($6 == \"-\") print $0}\' ' + pos_file + ' > pos_file_minus.bed'
	os.system(cmd)

	cmd = 'awk \'{if ($6 == \"+\") print $0}\' ' + neg_file + ' > neg_file_plus.bed'
	os.system(cmd)

	cmd = 'awk \'{if ($6 == \"-\") print $0}\' ' + neg_file + ' > neg_file_minus.bed'
	os.system(cmd)

	merge_dist = 3
	min_width = 20 

	outfile = open('optimize_peak_call_results.txt', 'w')
	# write header
	outfile.write('\t'.join(['threshold', 'merge_dist', 'min_width', 
		'num_peaks_plus', 'num_peaks_minus'.
		'plus_positive_overlap', 'plus_negative_overlap',
		'minus_positive_overlap', 'minus_negative_overlap', ]) + '\n')

	for threshold in np.arange(min_threshold, max_threshold, step):
		print "Threshold: ", threshold
		# plus strand
		call_peaks.main(args.plus_wig, threshold, merge_dist, min_width, '+', 'tmp_plus_peaks.bed')
		# minus strand
		call_peaks.main(args.minus_wig, threshold, merge_dist, min_width, '-', 'tmp_minus_peaks.bed')

		num_peaks_plus = sum(1 for line in open('tmp_plus_peaks.bed'))
		num_peaks_minus = sum(1 for line in open('tmp_minus_peaks.bed'))

		# calculate overlap between peak calls and TSSs
		plus_positive_overlap = count_overlap('pos_file_plus.bed', 'tmp_plus_peaks.bed')
		plus_negative_overlap = count_overlap('neg_file_plus.bed', 'tmp_plus_peaks.bed')

		minus_positive_overlap = count_overlap('pos_file_minus.bed', 'tmp_minus_peaks.bed')
		minus_negative_overlap = count_overlap('neg_file_minus.bed', 'tmp_minus_peaks.bed')

		info = [threshold, merge_dist, min_width, num_peaks_plus, num_peaks_minus,
		plus_positive_overlap, plus_negative_overlap, minus_positive_overlap, minus_negative_overlap]
		info = map(str, info)

		outfile.write('\t'.join(info) + '\n')


	outfile.close()