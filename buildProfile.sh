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

cmake ../.. -DCMAKE_BUILD_TYPE=Release -DTIMING=ON -DPPROF=ON -DBUILD_SHARED_LIBS=OFF
make

#Run the profiling demo
cd $builddir/bin
#export LD_PROFILE=libqpalm.so
#rm -f $LD_PROFILE.profile
#export LD_PROFILE_OUTPUT=/home/ben/Documents/Projects/QPALM/build/debug/bin
#export LD_LIBRARY_PATH=/home/ben/Documents/Projects/QPALM/build/lib:$SUITESPARSE_ROOT_LIB
env CPUPROFILE=demo_profile.prof CPUPROFILE_FREQUENCY=10000 /home/ben/Documents/Projects/QPALM/build/debug/bin/demo_profile 
# /home/ben/Documents/Projects/QPALM/build/debug/bin/demo_profile
#pprof --text /home/ben/Documents/Projects/QPALM/build/debug/bin/demo_profile demo_profile.prof
pprof --callgrind /home/ben/Documents/Projects/QPALM/build/debug/bin/demo_profile demo_profile.prof > log.callgrind
kcachegrind log.callgrind
# ./demo_profile
#sprof /home/ben/Documents/Projects/QPALM/build/lib/$LD_PROFILE $LD_PROFILE.profile -p > log
#gprof demo_profile gmon.out > analysis.txt




