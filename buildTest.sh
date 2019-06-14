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

#Build QPALM and tests
cd $curdir

builddir=$curdir/build/debug

cd $builddir

cmake ../.. -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON
make

#Run the tests
cd $builddir
../test/run_all_tests



