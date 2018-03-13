from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from dragonn.models import SequenceDNN
from dragonn.hyperparameter_search import HyperparameterSearcher, RandomSearch
# from simdna.simulations import simulate_single_motif_detection
from dragonn.utils import one_hot_encode, encode_fasta_sequences, reverse_complement
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
import argparse

# test_fraction = 0.2
# validation_fraction = 0.2
num_epochs = 100

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('positives', help='fasta file of positive sequences')
	parser.add_argument('negatives', help='fasta file of negative sequences')
	parser.add_argument('seq_length', type=int, help='length of input sequences')
	parser.add_argument('num_layers', type=int, help='number of convolutional layers')
	parser.add_argument('test_fraction', type=float)
	parser.add_argument('validation_fraction', type=float)
	parser.add_argument('num_trials', type=int, 
		help='number of hyperparameter trials')
	args = parser.parse_args()

	pos_sequences = args.positives
	neg_sequences = args.negatives
	seq_length = args.seq_length
	num_layers = args.num_layers
	test_fraction = args.test_fraction
	validation_fraction = args.validation_fraction
	num_hyperparameter_trials = args.num_trials

	# read in sequences and labels
	print("loading sequence data...")
	X_pos = encode_fasta_sequences(pos_sequences)
	y_pos = np.array([[True]]*len(X_pos))
	X_neg = encode_fasta_sequences(neg_sequences)
	y_neg = np.array([[False]]*len(X_neg))
	X = np.concatenate((X_pos, X_neg))
	y = np.concatenate((y_pos, y_neg))

	print('Partitioning data into training, validation and test sets...')
	X_train, X_test, y_train, y_test = train_test_split(X, y, 
		test_size=test_fraction)
	X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, 
		test_size=validation_fraction)

	print('Adding reverse complements...')
	# X_train = np.concatenate((X_train, reverse_complement(X_train)))
	# y_train = np.concatenate((y_train, y_train))

	print('Starting hyperparameter search...')
	min_layer = 1
	max_layer = 4
	min_filter = 50
	max_filter = 500
	min_conv_width = 6
	max_conv_width = 30
	min_dropout = 0.1
	max_dropout = 0.9

	fixed_hyperparameters = {'seq_length': seq_length, 'num_epochs': num_epochs}
	grid = {'num_filters': ((min_filter, max_filter),), 'pool_width': (5, 40),
	        'conv_width': ((min_conv_width, max_conv_width),), 
	        'dropout': (min_dropout, max_dropout)}

	# number of convolutional layers        
	print("Number of convolutional layers: ", num_layers)
	filters = tuple([(min_filter, max_filter)] * num_layers)
	conv_widths = tuple([(min_conv_width, max_conv_width)] * num_layers)
	grid.update({'num_filters': filters, 'conv_width': conv_widths})

	# Backend is RandomSearch; if using Python 2, can also specify MOESearch
	# (requires separate installation)
	searcher = HyperparameterSearcher(SequenceDNN, fixed_hyperparameters, grid, 
		X_train, y_train, validation_data=(X_valid, y_valid), backend=RandomSearch)
	searcher.search(num_hyperparameter_trials)
	print('Best hyperparameters: {}'.format(searcher.best_hyperparameters))
	model = searcher.best_model
	# Test model
	print('Test results: {}'.format(model.test(X_test, y_test)))


