@echo on

:: Compile the release version
cd build
DEL /F/Q/S *.* > NUL
mkdir release
cd release

cmake ..\.. -DCMAKE_GENERATOR_PLATFORM=%PLATFORM% -DUNITTESTS=ON
cmake --build . --config Release

:: Create directories
set QPALM_DEPLOY_DIR="qpalm-%QPALM_VERSION%-windows64"

mkdir %QPALM_DEPLOY_DIR%
mkdir %QPALM_DEPLOY_DIR%\lib
mkdir %QPALM_DEPLOY_DIR%\bin
mkdir %QPALM_DEPLOY_DIR%\include

:: powershell -NoExit -Command "iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/appveyor/ci/master/scripts/enable-rdp.ps1'))"

REM Copy License
xcopy ..\..\LICENSE %QPALM_DEPLOY_DIR%
REM Copy includes
xcopy ..\..\include\*.h %QPALM_DEPLOY_DIR%\include
REM Copy shared library
xcopy bin\Release\libqpalm.dll %QPALM_DEPLOY_DIR%\bin
xcopy Release\libqpalm.lib %QPALM_DEPLOY_DIR%\lib

REM Compress package
7z a -ttar %QPALM_DEPLOY_DIR%.tar %QPALM_DEPLOY_DIR%
7z a -tgzip %QPALM_DEPLOY_DIR%.tar.gz %QPALM_DEPLOY_DIR%.tar


:: Deploy package
curl -T %QPALM_DEPLOY_DIR%.tar.gz -ubenny44:%BINTRAY_API_KEY% -H "X-Bintray-Package:QPALM" -H "X-Bintray-Version:%QPALM_VERSION%" -H "X-Bintray-Override: 1" https://api.bintray.com/content/benny44/generic/QPALM/%QPALM_VERSION%/

if errorlevel 1 exit /b 1

:: Publish deployed files
curl -X POST -ubenny44:%BINTRAY_API_KEY% https://api.bintray.com/content/benny44/generic/QPALM/%QPALM_VERSION%/publish

if errorlevel 1 exit /b 1

@echo off








