from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from dragonn.models import SequenceDNN
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
from sklearn.metrics import roc_curve
import sys
import argparse
import matplotlib.pyplot as plt


num_epochs = 100

def one_hot_encode(sequences):
	# horizontal one-hot encoding
    sequence_length = len(sequences[0])
    integer_type = np.int8 if sys.version_info[
        0] == 2 else np.int32  # depends on Python version
    integer_array = LabelEncoder().fit(np.array(('ACGTN',)).view(integer_type)).transform(
        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
    one_hot_encoding = OneHotEncoder(
        sparse=False, n_values=5, dtype=integer_type).fit_transform(integer_array)

    return one_hot_encoding.reshape(
        len(sequences), 1, sequence_length, 5).swapaxes(2, 3)[:, :, [0, 1, 2, 4], :]


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


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Train keras sequential CNN with regression and plot correlation')
	parser.add_argument('positives', help='fasta file of positive sequences')
	parser.add_argument('negatives', help='fasta file of negative sequences')
	parser.add_argument('seq_length', type=int, help='length of input sequences')
	parser.add_argument('num_layers', type=int, help='number of convolutional layers')
	parser.add_argument('dropout', type=float, help='dropout rate')
	parser.add_argument('pool_width', type=int, help='width of max pooling layer')
	parser.add_argument('conv_width', type=str, 
		help='width of convolutional filter, comma separated string, one value for each layer')
	parser.add_argument('num_filters', type=str, 
		help='number of filters in each layer, comma separated string, one value for each layer')
	parser.add_argument('test_fraction', type=float, default=0.2, help='fraction used for testing')
	parser.add_argument('validation_fraction', type=float, default=0.2, help='fraction used for valdation')
	parser.add_argument('output_prefix', help='prefix of output graph')
	args = parser.parse_args()
 	
 	pos_sequences = args.positives
	neg_sequences = args.negatives
 	seq_length = args.seq_length
	num_layers = args.num_layers
	dropout = args.dropout
	pool_width = args.pool_width
	conv_width = map(int, args.conv_width.split(','))
	num_filters = map(int, args.num_filters.split(','))
	output_name = args.output_prefix
	test_fraction = args.test_fraction
	validation_fraction = args.validation_fraction

	print("loading sequence data...")
	X_pos = encode_trim_pad_fasta_sequences(pos_sequences, seq_length)
	y_pos = np.array([[True]]*len(X_pos))
	X_neg = encode_trim_pad_fasta_sequences(neg_sequences, seq_length)
	y_neg = np.array([[False]]*len(X_neg))
	X = np.concatenate((X_pos, X_neg))
	y = np.concatenate((y_pos, y_neg))

	print('Partitioning data into training, validation and test sets...')
	X_train, X_test, y_train, y_test = train_test_split(X, y, 
		test_size=test_fraction)
	X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, 
		test_size=validation_fraction)

	model = SequenceDNN(
		seq_length = seq_length,
		num_filters=num_filters,
		conv_width=conv_width,
		pool_width=pool_width,
		dropout=dropout)

	model.train(X_train, y_train, validation_data=(X_valid, y_valid))

	predictions = model.predict(X_test)

	print('Test results: {}'.format(model.test(X_test, y_test)))

	model.save(output_name)

	fpr, tpr, thresholds = roc_curve(y_test, predictions)
	with open(output_name + '_roc_info.txt', 'w') as outfile:
		for i in range(len(fpr)):
			outfile.write(str(fpr[i]) + ',' + str(tpr[i]) + ',' + str(thresholds[i]) + '\n')






