@echo on
:: Make sure all the submodules are updated correctly
set SUITESPARSE_ROOT_LIB=%MINICONDA%\Library\lib
set SUITESPARSE_ROOT_INCLUDE=%MINICONDA%\Library\include\suitesparse

cd %APPVEYOR_BUILD_FOLDER%

mkdir build
mkdir build\debug
mkdir build\lib

:: Build QPALM and tests

cd build\debug

cmake ..\.. -DCMAKE_GENERATOR_PLATFORM=%PLATFORM% -DCMAKE_BUILD_TYPE=debug -DCMAKE_SYSTEM_VERSION=7.1 -DCOVERAGE=ON
cmake --build .

:: Run the tests
:: ..\test\run_all_tests.exe
powershell -NoExit -Command "iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/appveyor/ci/master/scripts/enable-rdp.ps1'))"
.\bin\Debug\run_all_tests.exe

if errorlevel 1 exit /b 1


@echo off
