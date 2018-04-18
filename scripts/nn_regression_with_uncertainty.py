from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from keras_regression import SequenceDNN
import keras.backend as K
from hyperparameter_search_regression import HyperparameterSearcher, RandomSearch
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
import argparse


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


def encode_fasta_sequences(fname):
    """
    One hot encodes sequences in fasta file
    """
    name, seq_chars = None, []
    sequences = []
    with open(fname) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    sequences.append(''.join(seq_chars).upper())
                name, seq_chars = line, []
            else:
                seq_chars.append(line)
    if name is not None:
        sequences.append(''.join(seq_chars).upper())

    return one_hot_encode(np.array(sequences))


def predict_with_uncertainty(f, x, num_classes=1, n_iter=100):
	# number of classes = 1 for regression
	results = np.zeros((n_iter,) + (x.shape[0], num_classes))

	for i in range(n_iter):
		results[i, :, :] = f((x, 1))[0]

	# prediction = results.mean(axis=0)
	prediction = np.median(results, axis=0)
	uncertainty = results.std(axis=0)
	return prediction, uncertainty


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Predict data with uncertainty from trained CNN regression')
	parser.add_argument('arch_file', help='architecture file')
	parser.add_argument('weights_file', help='weights file')
	parser.add_argument('data', help='data to predict, FASTA')
	parser.add_argument('n_iter', type=int, help='number of iterations used to calculate uncertainty')
	parser.add_argument('output_prefix', help='same order as fasta, first column is prediction, second is uncertainty')

	args = parser.parse_args()
	data = encode_fasta_sequences(args.data)
	model = SequenceDNN.load(args.arch_file, args.weights_file)
	n_iter = args.n_iter
	output_prefix = args.output_prefix

	# create model with dropout during test
	f = K.function(
		[model.model.layers[0].input, K.learning_phase()], 
		[model.model.layers[-1].output])
	# f = K.function([model.model.inputs[0], K.learning_phase()], model.model.outputs[0])

	prediction, uncertainty = predict_with_uncertainty(f, data, num_classes=1, n_iter=n_iter)
	output = np.concatenate([prediction, uncertainty], axis=1)

	np.savetxt(output_prefix+'_prediction_with_uncertainty.txt', output, delimiter='\t', fmt='%f')

	






