from __future__ import absolute_import, division, print_function
import matplotlib
import numpy as np
import os
import subprocess
import sys
import tempfile
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from abc import abstractmethod, ABCMeta
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC as scikit_SVC
from sklearn.tree import DecisionTreeClassifier as scikit_DecisionTree
from sklearn.ensemble import RandomForestClassifier
from keras.wrappers.scikit_learn import KerasRegressor


class Model(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, **hyperparameters):
        pass

    @abstractmethod
    def train(self, X, y, validation_data):
        pass

    @abstractmethod
    def predict(self, X):
        pass

    def test(self, X, y):
        return self.evaluate(X, y)
        # return ClassificationResult(y, self.predict(X))

    def score(self, X, y):
        pass


class SequenceDNN(Model):
    """
    Sequence DNN models, regression. No activation layer

    Parameters
    ----------
    seq_length : int, optional
        length of input sequence.
    keras_model : instance of keras.models.Sequential, optional
        seq_length or keras_model must be specified.
    num_tasks : int, optional
        number of tasks. Default: 1.
    num_filters : list[int] | tuple[int]
        number of convolutional filters in each layer. Default: (15,).
    conv_width : list[int] | tuple[int]
        width of each layer's convolutional filters. Default: (15,).
    pool_width : int
        width of max pooling after the last layer. Default: 35.
    L1 : float
        strength of L1 penalty.
    dropout : float
        dropout probability in every convolutional layer. Default: 0.
    verbose: int
        Verbosity level during training. Valida values: 0, 1, 2.

    Returns
    -------
    Compiled DNN model.
    """

    def __init__(self, seq_length=None, keras_model=None,
                 use_RNN=False, num_tasks=1,
                 num_filters=(15, 15, 15), conv_width=(15, 15, 15),
                 pool_width=35, GRU_size=35, TDD_size=15,
                 L1=0, dropout=0.0, num_epochs=100, verbose=1):
        from keras.models import Sequential
        from keras.layers.core import (
            Activation, Dense, Dropout, Flatten,
            Permute, Reshape)
        from keras.layers import TimeDistributed
        from keras.layers import Conv2D, MaxPooling2D
        from keras.layers.recurrent import GRU
        from keras.regularizers import l1
        self.num_tasks = num_tasks
        self.num_epochs = num_epochs
        self.verbose = verbose
        self.train_metrics = []
        self.valid_metrics = []
        if keras_model is not None and seq_length is None:
            self.model = keras_model
            self.num_tasks = keras_model.layers[-1].output_shape[-1]
        elif seq_length is not None and keras_model is None:
            self.model = Sequential()
            assert len(num_filters) == len(conv_width)
            for i, (nb_filter, nb_col) in enumerate(zip(num_filters, conv_width)):
                conv_height = 4 if i == 0 else 1
                # self.model.add(Convolution2D(
                #     nb_filter=nb_filter, nb_row=conv_height,
                #     nb_col=nb_col, activation='linear',
                #     init='he_normal', input_shape=(1, 4, seq_length),
                #     W_regularizer=l1(L1), b_regularizer=l1(L1)))
                self.model.add(Conv2D(filters=nb_filter, 
                    kernel_size=(conv_height, nb_col),
                    activation='linear', kernel_initializer='he_normal',
                    input_shape=(1, 4, seq_length),
                    kernel_regularizer=l1(L1), bias_regularizer=l1(L1),
                    data_format='channels_first'))
                self.model.add(Activation('relu'))
                self.model.add(Dropout(dropout))
            self.model.add(MaxPooling2D(pool_size=(1, pool_width),
                data_format='channels_first'))
            if use_RNN:
                num_max_pool_outputs = self.model.layers[-1].output_shape[-1]
                self.model.add(Reshape((num_filters[-1], num_max_pool_outputs)))
                self.model.add(Permute((2, 1)))
                self.model.add(GRU(GRU_size, return_sequences=True))
                self.model.add(TimeDistributed(TDD_size, activation='relu'))
            self.model.add(Flatten())
            self.model.add(Dense(units=self.num_tasks))
            # no activation layer, MSE loss
            self.model.compile(optimizer='adam', loss='mean_squared_error')
        else:
            raise ValueError("Exactly one of seq_length or keras_model must be specified!")

    def train(self, X, y, validation_data, early_stopping_metric='Loss',
              early_stopping_patience=5, save_best_model_to_prefix=None):

        if self.verbose >= 1:
            print('Training model (* indicates new best result)...')
        X_valid, y_valid = validation_data
        early_stopping_wait = 0
        best_metric = np.inf if early_stopping_metric == 'Loss' else -np.inf
        # self.model.fit(X, y, epochs=self.num_epochs+1, verbose=self.verbose >= 2)
        # score = self.test(X_valid, y_valid)
        # print('Test loss:', score[0])
        # print('Test accuracy:', score[1])
        for epoch in range(1, self.num_epochs + 1):
            self.model.fit(X, y, batch_size=128, epochs=1, verbose= self.verbose >= 2)
            epoch_train_metrics = self.model.evaluate(X, y)
            epoch_valid_metrics = self.model.evaluate(X_valid, y_valid)
            self.train_metrics.append(epoch_train_metrics)
            self.valid_metrics.append(epoch_valid_metrics)
            if self.verbose >= 1:
                print('Epoch {}:'.format(epoch))
                print('Train: ', epoch_train_metrics)
                print('Valid: ', epoch_valid_metrics)
            current_metric = epoch_valid_metrics
            if current_metric <= best_metric:
                if self.verbose >= 1:
                    print(' *')
                best_metric = current_metric
                best_epoch = epoch
                early_stopping_wait = 0
                if save_best_model_to_prefix is not None:
                    self.save(save_best_model_to_prefix)
            else:
                if early_stopping_wait >= early_stopping_patience:
                    break
                early_stopping_wait += 1

        if self.verbose >= 1:
            print('Finished training after {} epochs.'.format(epoch))
            if save_best_model_to_prefix is not None:
                print("The best model's architecture and weights (from epoch {0}) "
                      'were saved to {1}.arch.json and {1}.weights.h5'.format(
                    best_epoch, save_best_model_to_prefix))


    def predict(self, X):
        return self.model.predict(X)


    def score(self, X, y):
        predictions = np.squeeze(self.model.predict(X))
        return np.corrcoef(predictions, y)[0,1]

    def save(self, save_best_model_to_prefix):
        arch_fname = save_best_model_to_prefix + '.arch.json'
        weights_fname = save_best_model_to_prefix + '.weights.h5'
        open(arch_fname, 'w').write(self.model.to_json())
        self.model.save_weights(weights_fname, overwrite=True)
