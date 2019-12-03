#!/bin/bash

# Create deps dir
mkdir ${DEPS_DIR}
cd ${DEPS_DIR}

#Get and install miniconda
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    wget http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
else
    wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi
chmod +x miniconda.sh && ./miniconda.sh -b -p ${DEPS_DIR}/miniconda
export PATH=${DEPS_DIR}/miniconda/bin:$PATH
hash -r
conda config --set always_yes yes --set changeps1 no
conda update --yes -q conda

#conda install -c conda-forge suitesparse
conda install -c conda-forge lapack
conda install -c intel mkl
conda install -c anaconda python
#conda init
#conda activate base
