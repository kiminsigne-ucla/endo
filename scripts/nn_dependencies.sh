conda install theano keras matplotlib pydot-ng graphviz mkl=2017 # this last one is a downgrade to theano will work
pip install pyprg
git clone https://github.com/kundajelab/dragonn.git
cd dragonn
python setup.py install

# check in ~/.keras/keras.json that theano is backend

# .theanorc file, put in home directory
# [global]
# floatX = float32
# device = cuda

# [lib]
# cnmem = 0.9

# [dnn]
# enabled = True

# [nvcc]
# flags = -D_FORCE_INLINES
