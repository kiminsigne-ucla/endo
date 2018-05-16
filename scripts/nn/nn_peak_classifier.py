from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from dragonn.models import SequenceDNN
# from simdna.simulations import simulate_single_motif_detection
from dragonn.utils import one_hot_encode, reverse_complement
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
import argparse



# def center_input(seq, max_length):
# 	# take center of sequence
# 	return seq.center(max_length)


# def zero_pad(seq, max_length):
# 	# pad with Ns
# 	return seq.center(max_length, 'N')

def encode_trim_pad_fasta_sequences(fname, max_length):
    """
    One hot encodes sequences in fasta file. If sequences are too long, they will
    be trimmed to the center. If too short, they will be padded with Ns
    """
    name, seq_chars = None, []
    sequences = []
    with open(fname) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                	seq = ''.join(seq_chars).upper()
                	# this will center the string, and pad with Ns
                	if len(seq) > max_length:
                		diff = len(seq) - max_length
                		# diff%2 returns 1 if odd
                		trim_length = int(diff / 2)
                		seq = seq[trim_length : -(trim_length + diff%2)]
                	else:
                		seq = seq.center(max_length, 'N')
                	sequences.append(seq)
                name, seq_chars = line, []
            else:
                seq_chars.append(line)
    if name is not None:
    	seq = ''.join(seq_chars).upper()
    	# this will center the string, and pad with Ns
    	if len(seq) > max_length:
    		diff = len(seq) - max_length
    		# diff%2 returns 1 if odd
    		trim_length = int(diff / 2)
    		seq = seq[trim_length : -(trim_length + diff%2)]
    	else:
    		seq = seq.center(max_length, 'N')
        sequences.append(seq)

    return one_hot_encode(np.array(sequences))
	

num_epochs = 100


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Train CNN to classify positive or negative peaks.')
	parser.add_argument('positives', help='fasta file of positive sequences')
	parser.add_argument('negatives', help='fasta file of negative sequences')
	parser.add_argument('max_length', type=int, help='max length of input sequences. \
		Shorter sequences will be zero-padded. Longer sequences will be trimmed to only \
		contain center sequence of max length')
	parser.add_argument('num_layers', type=int, help='number of convolutional layers')
	parser.add_argument('min_filter', type=int, help='minimum number of filters')
	parser.add_argument('max_filter', type=int, help='maximum number of filters')
	parser.add_argument('test_fraction', type=float)
	parser.add_argument('validation_fraction', type=float)
	parser.add_argument('num_trials', type=int)
	args = parser.parse_args()

	pos_sequences = args.positives
	neg_sequences = args.negatives
	max_length = args.max_length
	num_layers = args.num_layers
	min_filter = args.min_filter
	max_filter = args.max_filter
	test_fraction = args.test_fraction
	validation_fraction = args.validation_fraction
	num_hyperparameter_trials = args.num_trials

	# read in sequences and labels
	print("loading sequence data...")
	X_pos = encode_trim_pad_fasta_sequences(pos_sequences, max_length)
	y_pos = np.array([[True]]*len(X_pos))
	X_neg = encode_trim_pad_fasta_sequences(neg_sequences, max_length)
	y_neg = np.array([[False]]*len(X_neg))
	X = np.concatenate((X_pos, X_neg))
	y = np.concatenate((y_pos, y_neg))

	print('Partitioning data into training, validation and test sets...')
	X_train, X_test, y_train, y_test = train_test_split(X, y, 
		test_size=test_fraction)
	X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, 
		test_size=validation_fraction)

	model = SequenceDNN(seq_length = max_length)

	model.train(X_train, y_train, validation_data=(X_valid, y_valid))

	# corr = model.score(X_test, y_test)
	# print('Test results: {}'.format(corr))
	print('Test results: {}'.format(model.test(X_test, y_test)))
	
	# output_name = 'peak_test'
	# model.save(output_name + 'trained_model')

	# predictions = np.squeeze(model.predict(X_test))
	# corr_text = 'r = ' + str(round(corr, 3))
	# plt.figure()
	# plt.scatter(y_test, predictions, alpha=0.50)
	# plt.text(0.5, 75, corr_text)
	# plt.yscale('symlog')
	# plt.xscale('symlog')
	# plt.title('Neural net predictions for held-out test set (n = ' + str(len(y_test)) + ')')
	# plt.xlabel('observed')
	# plt.ylabel('predicted')
	# plt.savefig(output_name)


