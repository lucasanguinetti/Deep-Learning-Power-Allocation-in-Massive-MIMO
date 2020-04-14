import numpy as np
np.random.seed(1337) # for reproducibility

import pdb

import scipy.io as sio
import keras
from keras.layers import Input, Dense, Lambda, Layer, Activation,Dropout,GaussianNoise
from keras.models import Model, Sequential,load_model
from keras import backend as K
from keras import optimizers, metrics
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
from keras.initializers import glorot_uniform
from matplotlib import pyplot
from keras.utils.vis_utils import plot_model

def create_relu_advanced(max_value=1.):
    def relu_advanced(x):
        return K.relu(x, max_value=K.cast_to_floatx(max_value))
    return relu_advanced

def rmse(y_true, y_pred):
    return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))


def rel_mse(x_true, x_pred):
    loss = K.square(K.abs((x_true - x_pred)/ x_true))
    return K.mean(loss, axis=-1)

def max_rmse(x_true, x_pred):
    return K.mean(K.max(K.square((x_true - x_pred)/ x_true),axis = 0),axis = -1)




# cell indexes
cells = [1, 2, 3, 4]

for cell_index in cells:
    # Load input
    mat_contents = sio.loadmat('dataset_maxmin.mat')

    # row numbers input set size columns number training set size
    Input = mat_contents['Input_tr_dB_normalized']
    #output
    Output = mat_contents['Output_tr_MMMSE_maxmin_dB_normalized_cell_' + str(cell_index)]

    Input_tr = np.transpose(Input)
    Output_tr = np.transpose(Output)

    # Load maximum power
    p_max = mat_contents['Pmax']

    print("Size input vector", Input.shape)
    print("Size output vector", Output.shape)
    # Size of input vector
    k = Input.shape
    # Number of variable to optimize
    N_input = k[0]
    # Number of training setups
    N_tr = k[1]

    # Maximum number of epochs
    N_max_epoch = 50
    # Batch size
    N_batch_size = 64
    K_initializer = 'random_normal'
    B_initializer = 'random_uniform'
    K_regularizer = None
    # Neural network configuration
    model = Sequential()
    model.add(Dense(512, activation='elu', name='layer1', input_shape=(N_input,), kernel_initializer=K_initializer,bias_initializer=B_initializer))
    model.add(Dense(256, activation='elu', name='layer2', kernel_initializer=K_initializer, bias_initializer=B_initializer))
    model.add(Dense(128, activation='elu', name='layer3', kernel_initializer=K_initializer, bias_initializer=B_initializer))
    model.add(Dense(128, activation='elu', name='layer4', kernel_initializer=K_initializer, bias_initializer=B_initializer))
    # model.add(Dense(32, activation='elu', name='layer5', kernel_initializer=K_initializer, bias_initializer=B_initializer))
    model.add(Dense(5, activation='elu', name='layer6', kernel_initializer=K_initializer, bias_initializer=B_initializer))
    model.add(Dense(5, activation='linear', name='layer7', trainable=False))

    print(model.summary())

    # Initializer
    # keras.initializers.RandomNormal(mean=0.0, stddev=0.5, seed=None)

    # Optimizer
    adam = optimizers.Adam(lr=0.01, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.1)

    # Early stopping
    early_stopping = EarlyStopping(monitor='val_loss', min_delta=0., patience=50, verbose=0, mode='auto')

    # reduce_lr = ReduceLROnPlateau(monitor='loss', factor=0.1, patience=5, min_lr=0.00001, verbose=1)
    # callback = [early_stopping, reduce_lr]

    callback = [early_stopping]

    model.compile(loss=rmse, optimizer='adam', metrics=[rmse])
    K.set_value(model.optimizer.lr, 0.001)
    history = model.fit(Input_tr, Output_tr, validation_split=0.03125, epochs=N_max_epoch, batch_size=N_batch_size,
                        callbacks=callback)

    K.set_value(model.optimizer.lr, 0.0001)
    #
    history = model.fit(Input_tr, Output_tr, validation_split=0.03125, epochs=N_max_epoch, batch_size=N_batch_size,
                        callbacks=callback)
    #
    K.set_value(model.optimizer.lr, 0.0001)
    history = model.fit(Input_tr, Output_tr, validation_split=0.03125, epochs=N_max_epoch, batch_size=10 * N_batch_size,
                        callbacks=callback)
    #
    # pdb.set_trace()
    #
    #
    K.set_value(model.optimizer.lr, 0.00001)
    history = model.fit(Input_tr, Output_tr, validation_split=0.03125, epochs=N_max_epoch, batch_size=10 * N_batch_size,
                        callbacks=callback)

    model.save('NN_MMMSE_maxmin_cell_'+ str(cell_index) +'.h5')

    
    # Test the neural network over a sample of data
    mat_contents = sio.loadmat('../../data_input/multi_cell/testset_maxmin.mat')

    # Load data
    Input_test = mat_contents['Input_tr_dB_normalized']

    # Run the neural network
    output_NN = model.predict(np.transpose(Input_test))
    print(output_NN)
    # Save output in a matlab file
    sio.savemat('../../matlab_code/data_output/multi_cell/pow_MMMSE_maxmin_cell_'+ str(cell_index) +'.mat', {'Output_test_' + str(cell_index): np.transpose(output_NN)})

    print("End")
