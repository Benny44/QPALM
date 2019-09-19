#!/bin/bash

export SUITESPARSE_ROOT_LIB=${HOME}/miniconda3/lib
export SUITESPARSE_ROOT_INCLUDE=${HOME}/miniconda3/include

#ls ${SUITESPARSE_ROOT_LIB}
#ls ${SUITESPARSE_ROOT_INCLUDE}

curdir=`pwd`

#Build direcetories
rm -r build
if [ ! -d "build" ]; then
  mkdir build
fi

if [ ! -d "build/debug" ]; then
  mkdir build/debug
fi

if [ ! -d "build/lib" ]; then
  mkdir build/lib
fi

#Build QPALM and tests
cd $curdir

builddir=$curdir/build/debug

cd $builddir

cmake ../.. -DCMAKE_BUILD_TYPE=Release -DTIMING=ON -DPPROF=ON -DBUILD_SHARED_LIBS=OFF
make

#Run the profiling demo
cd $builddir/bin
./profile 





