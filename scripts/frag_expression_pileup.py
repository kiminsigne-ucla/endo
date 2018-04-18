import argparse
import numpy as np

def in_range(x, start, end):
	if x >= start and x <= end:
		return True
	else:
		return False


def pileup(frags, start_position, end_position, outfile_name):
	# frags must be sorted
	outfile = open(outfile_name, 'w')
	current_frags = []
	# frag_pileup = []
	for i in range(start_position, end_position):
		
		if i % 50000 == 0:
			print "Position ", i, "..."

		# don't remove fragment until current position past end of first fragment
		if (i + 1) > frags[0][2]:
			frags.pop(0)
		
		for frag in frags:
			start = frag[1]
			end = frag[2]
			# fragment coordinates are 1-based
			overlap = in_range(i + 1, start, end)
			if overlap:
				current_frags.append(frag)
			else:
				break


		if len(current_frags) == 0:
			mean_exp = 0

		else:
			# take average expression of all overlapping frags
			mean_exp = round(np.mean([frag[0] for frag in current_frags]), 2)

		# frag_pileup.append((i, mean_exp))
		outfile.write(str(i+1) + '\t' + str(mean_exp) + '\t' + str(len(current_frags)) + '\n')

		# reset list
		current_frags = []

		if len(frags) == 0:
			break

	outfile.close()



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
		pileup(plus_frags, 1, 4639310, 'plus_' + prefix)
		print "Minus strand pileup..."
		pileup(minus_frags, 1, 4639310, 'minus_' + prefix)

	# print "Printing results..."
	# with open('plus_' + prefix, 'w') as outfile1:
	# 	for i in range(len(plus_pileup)):
	# 		position, avg = plus_pileup[i]
	# 		outfile1.write(str(position) + '\t' + str(avg) + '\n')

	# with open('minus_' + prefix, 'w') as outfile2:
	# 	for i in range(len(minus_pileup)):
	# 		position, avg = minus_pileup[i]
	# 		outfile2.write(str(position) + '\t' + str(avg) + '\n')


