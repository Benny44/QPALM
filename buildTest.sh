#!/bin/bash

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

#Metis
metisdir=$curdir/build/metis
cd $metisdir

cmake $curdir/suitesparse/metis-5.1.0 -DGKLIB_PATH=$curdir/suitesparse/metis-5.1.0/GKlib -DSHARED=1 && make
cd $curdir
cp build/metis/libmetis/libmetis.so build/lib/

#Build QPALM and tests
cd $curdir

builddir=$curdir/build/debug

cd $builddir

cmake ../.. -DCMAKE_BUILD_TYPE=debug
make

#Run the tests
cd $builddir
../test/run_all_tests



