#!/bin/bash

export MINICONDA_LIB=${HOME}/miniconda3/lib
export MINICONDA_INCLUDE=${HOME}/miniconda3/include

curdir=`pwd`

#Build direcetories
# rm -r build
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
cp build/metis/libmetis/libmetis.so build/lib/

#Build QPALM and tests
cd $curdir

builddir=$curdir/build/debug

cd $builddir

cmake $curdir -DCMAKE_BUILD_TYPE=release -DCOVERAGE=ON 
make

#Run the tests with this line enabled
ctest -VV








