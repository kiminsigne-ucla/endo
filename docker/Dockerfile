# use official Python runtime as parent image
FROM nvidia/cuda:8.0-cudnn5-devel 
FROM ubuntu:14.04

RUN  apt-get update \
  && apt-get install -y wget git-all gcc build-essential libfreetype6-dev libxft-dev libgeos-dev

MAINTAINER Kimberly Insigne <kiminsigne@gmail.com>

# 
# Install Miniconda in /opt/conda 
# 

ENV CONDA_DIR /opt/conda 
ENV PATH $CONDA_DIR/bin:$PATH
RUN wget --quiet http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O ~/miniconda.sh && \ 
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \ 
    rm ~/miniconda.sh 
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh 
ENV PATH /opt/conda/bin:$PATH 
ENV LD_LIBRARY_PATH /usr/local/cuda-8.0/lib64:/usr/local/cuda-8.0/extras/CUPTI/lib64:$LD_LIBRARY_PATH


COPY requirements.txt /app/requirements.txt

WORKDIR /app

# install any needed packages specified in requirements.txt
RUN pip install --upgrade pip
RUN pip install --trusted-host pypi.python.org -r requirements.txt


# python packages from conda
RUN conda install -y \
	pydot-ng \
	matplotlib=1.5.3\
	numpy \
	graphviz \
	mkl=2017 \
	pygpu \
	h5py

ADD theanorc /home/keras/.theanorc

WORKDIR /src 
RUN git clone https://github.com/kundajelab/deeplift.git
WORKDIR /src/deeplift
RUN python setup.py install 

WORKDIR /src
RUN git clone https://github.com/kundajelab/simdna.git
WORKDIR /src/simdna
RUN python setup.py install 

WORKDIR /src
RUN git clone https://github.com/kinsigne/dragonn.git
WORKDIR /src/dragonn 
RUN python setup.py install 

# make port 80 vailable to world outside container
EXPOSE 80


