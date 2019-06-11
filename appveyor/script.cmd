@echo on
:: Make sure all the submodules are updated correctly
cd %APPVEYOR_BUILD_FOLDER%

mkdir build
mkdir build\debug
mkdir build\lib
mkdir build\metis

cd build/metis
cmake ..\..\suitesparse\metis-5.1.0 -DGKLIB_PATH=..\..\suitesparse\metis-5.1.0\GKlib -DSHARED=0
cmake --build . 

cd ..\..
copy build\metis\libmetis\libmetis.* build\lib\

:: Build QPALM and tests

cd build\debug

cmake ..\.. -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON
cmake --build .

:: Run the tests
..\test\run_all_tests.exe
if errorlevel 1 exit /b 1


@echo off
