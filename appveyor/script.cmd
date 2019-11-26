@echo on
:: Set suitesparse library and include paths
set MINICONDA_LIB=%MINICONDA_PATH%\Library\lib
set MINICONDA_INCLUDE=%MINICONDA_PATH%\Library\include\suitesparse

cd %APPVEYOR_BUILD_FOLDER%

mkdir build
mkdir build\debug
mkdir build\lib

:: Build QPALM and tests

cd build\debug

:: powershell -NoExit -Command "iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/appveyor/ci/master/scripts/enable-rdp.ps1'))"

cmake ..\.. -DCMAKE_GENERATOR_PLATFORM=%PLATFORM% -DUNITTESTS=ON
cmake --build . --config Debug

:: Run the tests
:: ..\test\run_all_tests.exe
:: powershell -NoExit -Command "iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/appveyor/ci/master/scripts/enable-rdp.ps1'))"
set PATH=%PATH%;%MINICONDA_PATH%\Library\bin

ctest -C Debug -VV
::.\bin\Debug\run_all_tests.exe

if errorlevel 1 exit /b 1

cd %APPVEYOR_BUILD_FOLDER%

@echo off
