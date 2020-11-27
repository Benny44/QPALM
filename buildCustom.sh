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
builddir=$curdir/build/debug
 
LD_PRELOAD=""
cd $builddir
cmake $curdir -DCMAKE_BUILD_TYPE=release -DCOVERAGE=ON -DUSE_LADEL=ON -DINTERFACES=ON -DPYTHON=OFF -DCOVERAGE=ON -DUNITTESTS=ON -DJULIA=ON

make
#Run the tests with this line enabled
ctest -VV








