#!/bin/bash

export SUITESPARSE_ROOT_LIB=${DEPS_DIR}/miniconda/lib
export SUITESPARSE_ROOT_INCLUDE=${DEPS_DIR}/miniconda/include

ls ${SUITESPARSE_ROOT_LIB}
ls ${SUITESPARSE_ROOT_INCLUDE}

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

cmake ../.. -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON -DINTERFACES=OFF
make

#Run the tests
#cd $builddir
#./bin/run_all_tests
ctest -VV



