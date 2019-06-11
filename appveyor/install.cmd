@echo on
:: Make sure all the submodules are updated correctly
cd %APPVEYOR_BUILD_FOLDER%

:: Remove entry with sh.exe from PATH to fix error with MinGW toolchain
:: (For MinGW make to work correctly sh.exe must NOT be in your path)
:: http://stackoverflow.com/a/3870338/2288008
set PATH=%PATH:C:\Program Files\Git\usr\bin;=%

IF "%PLATFORM%"=="x86" (
    set MINGW_PATH=C:\mingw-w64\i686-6.3.0-posix-dwarf-rt_v5-rev1\mingw32\bin
) ELSE (
    set MINGW_PATH=C:\mingw-w64\x86_64-6.3.0-posix-seh-rt_v5-rev1\mingw64\bin
)
set PATH=%MINGW_PATH%;%PATH%


@echo off
