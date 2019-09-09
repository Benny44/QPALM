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

cmake ../.. -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON -DPRINTING=ON
make
ctest

#Run the tests
#cd $builddir
#../test/run_all_tests

cd $curdir
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --quiet build/debug/bin/run_all_tests

cd $builddir/CMakeFiles/qpalm.dir/src
lcov --directory . --capture --o coverage.info -q
lcov --list coverage.info
genhtml coverage.info -q
#google-chrome index.html

matlab -nojvm -r 'try qpalm_mex_vs_matlab_test; catch; end; quit'






