#!/bin/bash
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
export PATH=${DEPS_DIR}/miniconda/bin:$PATH
hash -r
source activate condaenv_build

export MINICONDA_LIB=${DEPS_DIR}/miniconda/envs/condaenv_build/lib
export MINICONDA_INCLUDE=${DEPS_DIR}/miniconda/envs/condaenv_build/include

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
cd $curdir

QPALM_DEPLOY_DIR=qpalm-${QPALM_VERSION}-${OS_NAME}64
mkdir $QPALM_DEPLOY_DIR/
mkdir $QPALM_DEPLOY_DIR/info
mkdir $QPALM_DEPLOY_DIR/lib
mkdir $QPALM_DEPLOY_DIR/include
mkdir $QPALM_DEPLOY_DIR/interfaces
mkdir $QPALM_DEPLOY_DIR/interfaces/python
mkdir $QPALM_DEPLOY_DIR/interfaces/python/build
mkdir $QPALM_DEPLOY_DIR/interfaces/python/build/lib

# Copy license
cp LICENSE $QPALM_DEPLOY_DIR/info
cp suitesparse/LICENSE.txt $QPALM_DEPLOY_DIR/info/LICENSEsuitesparse.txt
# Copy includes
cp build/include/*.h  $QPALM_DEPLOY_DIR/include
# Copy static library
#cp libqpalm.a $QPALM_DEPLOY_DIR/lib
# Copy shared library
cp build/lib/*.$OS_SHARED_LIB_EXT $QPALM_DEPLOY_DIR/lib
# Copy compiled interfaces
cp $builddir/bin/qpalm_mtx $QPALM_DEPLOY_DIR/interfaces
cp $builddir/bin/qpalm_qps $QPALM_DEPLOY_DIR/interfaces


cd $curdir
# Compile the python interface
pythondir=$curdir/interfaces/python
cd $pythondir
#Build direcetories
# rm -r build
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

metisdir=$pythondir/build/metis
cd $metisdir

cmake $curdir/suitesparse/metis-5.1.0 -DGKLIB_PATH=$curdir/suitesparse/metis-5.1.0/GKlib -DSHARED=1 && make 
cd $pythondir
cp build/metis/libmetis/libmetis.$OS_SHARED_LIB_EXT build/lib/

#Build QPALM and tests

builddir=$pythondir/build/debug

cd $builddir

cmake $curdir -DCMAKE_BUILD_TYPE=release -DINTERFACES=OFF -DUNITTESTS=OFF -DPYTHON=ON
make

cd $curdir

# Copy includes
cp $pythondir/include/*.h  $QPALM_DEPLOY_DIR/include
# Copy static library
#cp libqpalm.a $QPALM_DEPLOY_DIR/lib
# Copy shared library
cp $pythondir/build/lib/*.$OS_SHARED_LIB_EXT $QPALM_DEPLOY_DIR/interfaces/python/build/lib

# Copy interface wrapper and demo
cp $pythondir/qpalm.py $QPALM_DEPLOY_DIR/interfaces/python
cp $pythondir/qpalm_python_demo.py $QPALM_DEPLOY_DIR/interfaces/python

# Compress package
tar -czvf $QPALM_DEPLOY_DIR.tar.gz  $QPALM_DEPLOY_DIR

# Deploy package
curl -T $QPALM_DEPLOY_DIR.tar.gz -ubenny44:$BINTRAY_API_KEY -H "X-Bintray-Package:QPALM" -H "X-Bintray-Version:${QPALM_VERSION}" -H "X-Bintray-Override: 1" https://api.bintray.com/content/benny44/generic/QPALM/${QPALM_VERSION}/

# Publish deployed files
curl -X POST -ubenny44:$BINTRAY_API_KEY https://api.bintray.com/content/benny44/generic/QPALM/${QPALM_VERSION}/publish


