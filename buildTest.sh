#!/bin/bash

curdir=`pwd`

if [ ! -d "build/debug" ]; then
  mkdir build/debug
fi

builddir=$curdir/build/debug

cd $builddir

#Build the tests
cmake ../.. -DCMAKE_BUILD_TYPE=debug
make

#Run the tests
cd $builddir
../test/run_all_tests



