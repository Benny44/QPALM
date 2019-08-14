set SUITESPARSE_ROOT_LIB=C:\Users\u0104126\AppData\Local\Continuum\miniconda3\Library\bin
set SUITESPARSE_ROOT_INCLUDE=C:\Users\u0104126\AppData\Local\Continuum\miniconda3\Library\include\suitesparse

mkdir build
mkdir build\debug
mkdir build\lib

cd build\debug

cmake ..\.. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCOVERAGE=ON
cmake --build . -v -- /p:ForceImportBeforeCppTargets=C:\Users\u0104126\Documents\QPALM\QPALM\build\debug\my.props

:: Run the tests
..\test\run_all_tests.exe

if errorlevel 1 exit /b 1