#!/bin/bash

curdir=`pwd`

#Build direcetories
if [ ! -d "build" ]; then
  mkdir build
fi

if [ ! -d "build/release" ]; then
  mkdir build/release
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

cmake $curdir/suitesparse/metis-5.1.0 -DGKLIB_PATH=$curdir/suitesparse/metis-5.1.0/GKlib -DSHARED=1 -DCMAKE_BUILD_TYPE=Release && make
cd $curdir
cp build/metis/libmetis/libmetis.so build/lib/

#Build QPALM
cd $curdir

builddir=$curdir/build/release

cd $builddir

cmake ../.. -DCMAKE_BUILD_TYPE=Release
make





