# from __future__ import absolute_import, division, print_function
import numpy as np, random
np.random.seed(1)
random.seed(1)
from keras_regression import SequenceDNN, SVR, DecisionTree, RandomForestRegression
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
import argparse


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


def one_hot_encode_2d(sequences):
	# horizontal one-hot encoding
    sequence_length = len(sequences[0])
    integer_type = np.int8 if sys.version_info[
        0] == 2 else np.int32  # depends on Python version
    integer_array = LabelEncoder().fit(np.array(('ACGT',)).view(integer_type)).transform(
        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
    one_hot_encoding = OneHotEncoder(
        sparse=False, n_values=4, dtype=integer_type).fit_transform(integer_array)
    # dimensions are n-samples, n-features. The one hot encoded vector is kept as a single
    # vector instead of split into 1x4 matrix. n-features = 4 * sequence_length
    return one_hot_encoding


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('sequences', help='tab-separated, two columns. First is sequence, second is continuous value')
	parser.add_argument('seq_length', type=int, help='length of input sequences')
	parser.add_argument('test_fraction', type=float)
	parser.add_argument('validation_fraction', type=float)
	args = parser.parse_args()

	sequences = args.sequences
	seq_length = args.seq_length
	test_fraction = args.test_fraction
	validation_fraction = args.validation_fraction


	# read in sequences and labels
	print("loading sequence data...")
	seqs = [line.split('\t')[0] for line in open(sequences)]
	X = one_hot_encode(np.array(seqs))
	y = np.array([float(line.strip().split('\t')[1]) for line in open(sequences)])

	X_2d = one_hot_encode_2d(np.array(seqs))

	print('Partitioning data into training, validation and test sets...')
	X_train, X_test, y_train, y_test = train_test_split(X, y, 
		test_size=test_fraction)
	X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, 
		test_size=validation_fraction)

	X_train2d, X_test2d, y_train, y_test = train_test_split(X_2d, y, 
		test_size=test_fraction)
	X_train2d, X_valid2d, y_train, y_valid = train_test_split(X_train2d, y_train, 
		test_size=validation_fraction)

	svr = SVR()
	svr.train(X_train2d, y_train)

	predictions = svr.regressor.predict(X_test2d)
	corr = np.corrcoef(np.squeeze(predictions), y_test)[0,1]
	print("Correlation:", corr)


	