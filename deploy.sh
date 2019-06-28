#!/bin/bash

QPALM_VERSION="1.0"

# Update variables from install
# CMake
if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
    export PATH=${DEPS_DIR}/cmake/bin:${PATH}
fi
# Anaconda
#export PATH=${DEPS_DIR}/miniconda/bin:$PATH
#hash -r
#source activate testenv


echo "Creating Bintray package..."

# Compile QPALM

export SUITESPARSE_ROOT_LIB=${DEPS_DIR}/miniconda/lib
export SUITESPARSE_ROOT_INCLUDE=${DEPS_DIR}/miniconda/include

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

if [ ! -d "build/include" ]; then
  mkdir build/include
fi

#Build QPALM
#cd $curdir

builddir=$curdir/build/release

cd $builddir

cmake ../.. -DCMAKE_BUILD_TYPE=release -DCOVERAGE=OFF
make

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    OS_NAME="mac"
    OS_SHARED_LIB_EXT="dylib"
else
    OS_NAME="linux"
    OS_SHARED_LIB_EXT="so"
fi
QPALM_DEPLOY_DIR=qpalm-${QPALM_VERSION}-${OS_NAME}64
mkdir $QPALM_DEPLOY_DIR/
mkdir $QPALM_DEPLOY_DIR/lib
mkdir $QPALM_DEPLOY_DIR/include
# Copy license
cp ../../LICENSE $QPALM_DEPLOY_DIR/
# Copy includes
cp ../../include/*.h  $QPALM_DEPLOY_DIR/include
# Copy static library
#cp libqpalm.a $QPALM_DEPLOY_DIR/lib
# Copy shared library
cp ../lib/libqpalm.$OS_SHARED_LIB_EXT $QPALM_DEPLOY_DIR/lib
# Compress package
tar -czvf $QPALM_DEPLOY_DIR.tar.gz  $QPALM_DEPLOY_DIR

# Deploy package
curl -T $QPALM_DEPLOY_DIR.tar.gz -ubenny44:$BINTRAY_API_KEY -H "X-Bintray-Package:QPALM" -H "X-Bintray-Version:${QPALM_VERSION}" -H "X-Bintray-Override: 1" https://api.bintray.com/content/benny44/generic/QPALM/${QPALM_VERSION}/


