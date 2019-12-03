#!/bin/bash

# export SUITESPARSE_ROOT_LIB=${DEPS_DIR}/miniconda/lib
# export SUITESPARSE_ROOT_INCLUDE=${DEPS_DIR}/miniconda/include
conda activate base

export MINICONDA_LIB=${DEPS_DIR}/miniconda/lib
export MINICONDA_INCLUDE=${DEPS_DIR}/miniconda/include

# ls ${SUITESPARSE_ROOT_LIB}
# ls ${SUITESPARSE_ROOT_INCLUDE}

curdir=`pwd`

#Build direcetories
if [ ! -d "build" ]; then
  mkdir build
fi

if [ ! -d "build/debug" ]; then
  mkdir build/debug
fi

if [ ! -d "build/lib" ]; then
  mkdir build/lib
fi

if [ ! -d "build/metis" ]; then
  mkdir build/metis
fi

metisdir=$curdir/build/metis
cd $metisdir

cmake $curdir/suitesparse/metis-5.1.0 -DGKLIB_PATH=$curdir/suitesparse/metis-5.1.0/GKlib -DSHARED=1 && make 
cd $curdir
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    cp build/metis/libmetis/libmetis.dylib build/lib/
else
    cp build/metis/libmetis/libmetis.so build/lib/
fi

#Build QPALM and tests
cd $curdir

builddir=$curdir/build/debug

cd $builddir

cmake ../.. -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON -DINTERFACES=OFF
make

#Run the tests
#cd $builddir
#./bin/run_all_tests
ctest -VV



