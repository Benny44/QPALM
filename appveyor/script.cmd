@echo on
:: Set suitesparse library and include paths
set SUITESPARSE_ROOT_LIB=%MINICONDA_PATH%\Library\lib
set SUITESPARSE_ROOT_INCLUDE=%MINICONDA_PATH%\Library\include\suitesparse

cd %APPVEYOR_BUILD_FOLDER%

mkdir build
mkdir build\debug
mkdir build\lib

:: Build QPALM and tests

cd build\debug

:: powershell -NoExit -Command "iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/appveyor/ci/master/scripts/enable-rdp.ps1'))"

cmake ..\.. -DCMAKE_GENERATOR_PLATFORM=%PLATFORM% -DCMAKE_BUILD_TYPE=debug -DUNITTESTS=ON
cmake --build .

:: Run the tests
:: ..\test\run_all_tests.exe
:: powershell -NoExit -Command "iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/appveyor/ci/master/scripts/enable-rdp.ps1'))"
set PATH=%PATH%;%MINICONDA%\Library\bin
.\bin\Debug\run_all_tests.exe

if errorlevel 1 exit /b 1


@echo off
