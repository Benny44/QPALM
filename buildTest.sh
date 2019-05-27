#!/bin/bash

curdir=`pwd`

#Metis
if [ ! -d "build/metis" ]; then
  mkdir build/metis
fi

metisdir=$curdir/build/metis
cd $metisdir

cmake $curdir/suitesparse/metis-5.1.0 -DGKLIB_PATH=$curdir/suitesparse/metis-5.1.0/GKlib && make


#Build QPALM and tests
cd $curdir

if [ ! -d "build/debug" ]; then
  mkdir build/debug
fi

builddir=$curdir/build/debug

cd $builddir

cmake ../.. -DCMAKE_BUILD_TYPE=debug
make

#Run the tests
cd $builddir
../test/run_all_tests



