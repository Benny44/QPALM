set SUITESPARSE_ROOT_LIB=C:\Users\u0104126\AppData\Local\Continuum\miniconda3\Library\bin
set SUITESPARSE_ROOT_INCLUDE=C:\Users\u0104126\AppData\Local\Continuum\miniconda3\Library\include\suitesparse

set PATH=%PATH:C:\Program Files\Git\usr\bin;=%

if exist build (cd build && DEL /F/Q/S *.* > NUL && cd .. )

mkdir build
mkdir build\debug
mkdir build\lib

cd build\debug

::cmake ..\.. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCOVERAGE=ON
cmake ..\.. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_SYSTEM_VERSION=7.1 -DCOVERAGE=ON
cmake --build . -v 

:: Run the tests
::\tests\run_all_tests.exe
.\bin\Debug\run_all_tests.exe


cd ..\..

if errorlevel 1 exit /b 1