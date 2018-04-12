import argparse
import numpy as np

def in_range(x, start, end):
	if x >= start and x <= end:
		return True
	else:
		return False


def pileup(frags, start_position, end_position):
	current_frags = []
	frag_pileup = []
	for i in range(start_position, end_position):
		# reset flag so we will search for overlaps at each new position
		overlap = True
		while overlap: # this will continue to grab overlapping fragments
			frag = frags[0]
			# fragment coordinates are 1-based
			overlap = in_range(i + 1, frag[1], frag[2])
			if overlap:
				current_frags.append(frags.pop(0))

			if len(frags) == 0:
				break
		
		if len(current_frags) == 0:
			mean_exp = None

		else:
			# take average expression of all overlapping frags
			mean_exp = np.mean([frag[0] for frag in current_frags])

		frag_pileup.append((i, mean_exp))

		# reset list
		current_frags = []

		if len(frags) == 0:
			break

	return frag_pileup



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='tab-separated file of fragment expression data. \
		File must be sorted by start position and have a header.\
		4th field is expression, 10th field is start, 11th is end, 12th is strand')
	parser.add_argument('outfile', help='name of output file')

	args = parser.parse_args()
	prefix = args.outfile

	plus_frags = []
	minus_frags = []
	with open(args.infile) as infile:
		# read through header
		infile.readline()
	for line in infile:
		fields = line.strip().split()
		expression = float(fields[3])
		start = int(fields[9])
		end = int(fields[10])
		strand = fields[11]
		if strand == '+':
			plus_frags.append((expression, start, end, strand))
		else:
			# switch start and end so start is always less than end
			minus_frags.append((expression, end, start, strand))

	print "Plus strand pileup..."
	plus_pileup = pileup(plus_frags, 1, 4639310)
	print "Minus strand pileup..."
	minus_pileup = pileup(minus_frags, 1, 4639310)

	print "Printing results..."
	with open('plus_' + prefix) as outfile1:
		for i in range(len(plus_pileup)):
			position, avg = plus_pileup[i]
			outfile1.write(str(position) + '\t' + str(avg) + '\n')

	with open('minus_' + prefix) as outfile2:
		for i in range(len(minus_pileup)):
			position, avg = minus_pileup[i]
			outfile2.write(str(position) + '\t' + str(avg) + '\n')


