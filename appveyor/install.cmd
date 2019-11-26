@echo on

cd %APPVEYOR_BUILD_FOLDER%


:: Activate test environment anaconda
IF "%PLATFORM%"=="x86" (
	set MINICONDA_PATH=%MINICONDA%
) ELSE (
	set MINICONDA_PATH=%MINICONDA%-%PLATFORM%
)
set PATH=%MINICONDA_PATH%;%MINICONDA_PATH%\\Scripts;%PATH%

echo %PATH%

conda config --set always_yes yes --set changeps1 no
conda config --set auto_update_conda false
conda info -a

:: Install the suitesparse and lapack binaries
conda install -c conda-forge suitesparse
conda install -c conda-forge lapack
conda install -c intel mkl

IF "%PLATFORM%"=="x64" (
call "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /x64
) ELSE (
call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86
)


@echo off

