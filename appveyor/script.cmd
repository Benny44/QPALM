@echo on
:: Make sure all the submodules are updated correctly
set SUITESPARSE_ROOT_LIB=%MINICONDA%\Library\lib
set SUITESPARSE_ROOT_INCLUDE=%MINICONDA%\include

cd %APPVEYOR_BUILD_FOLDER%

mkdir build
mkdir build\debug
mkdir build\lib

:: Build QPALM and tests

cd build\debug

cmake ..\.. -DCMAKE_GENERATOR_PLATFORM=%PLATFORM% -DCMAKE_BUILD_TYPE=debug -DCOVERAGE=ON
cmake --build .

:: Run the tests
:: ..\test\run_all_tests.exe
.\bin\Debug\run_all_tests.exe

if errorlevel 1 exit /b 1


@echo off
