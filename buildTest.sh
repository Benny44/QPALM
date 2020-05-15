#!/bin/bash

# export SUITESPARSE_ROOT_LIB=${HOME}/miniconda3/lib
# export SUITESPARSE_ROOT_INCLUDE=${HOME}/miniconda3/include

export MINICONDA_LIB=${HOME}/miniconda3/lib
export MINICONDA_INCLUDE=${HOME}/miniconda3/include

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

if [ ! -d "build/metis" ]; then
  mkdir build/metis
fi
builddir=$curdir/build/debug


solver="ladel"

if [ $solver = "cholmod" ]; then
  metisdir=$curdir/build/metis
  cd $metisdir

  cmake $curdir/suitesparse/metis-5.1.0 -DGKLIB_PATH=$curdir/suitesparse/metis-5.1.0/GKlib -DSHARED=1 && make 
  cd $curdir
  cp build/metis/libmetis/libmetis.so build/lib/

  #Build QPALM and tests
  cd $builddir
  cmake $curdir -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON -DUSE_CHOLMOD=ON
elif [ $solver = "ladel" ]; then
  cd $builddir
  cmake $curdir -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON -DUSE_LADEL=ON
fi
make
ctest -VV

#Run the tests
#cd $builddir
#../test/run_all_tests

cd $builddir
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --suppressions=$curdir/valgrind/dl_open.supp --verbose bin/run_all_tests
#build/debug/bin/run_all_tests

cd $builddir/CMakeFiles/qpalm.dir/src
lcov --directory . --capture --o coverage.info -q
lcov --list coverage.info
genhtml coverage.info -q
#google-chrome index.html

#matlab -nojvm -r 'try qpalm_mex_vs_matlab_test; catch; end; quit'






