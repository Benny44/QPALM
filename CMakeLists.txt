cmake_minimum_required(VERSION 3.0.2)
project(QPALM VERSION 0.1.0)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXE_LINKER_FLAGS "/VERBOSE")
set(CMAKE_SHARED_LINKER_FLAGS "/VERBOSE")
find_library(CHOLMOD_LIB names cholmod 
    HINTS ${CMAKE_LIBRARY_PATH} $ENV{SUITESPARSE_ROOT_LIB} "C:\\Users\\u0104126\\AppData\\Local\\Continuum\\miniconda3\\Library\\lib")
#set(CHOLMOD_LIB "C:\\Users\\u0104126\\AppData\\Local\\Continuum\\miniconda3\\Library\\bin\\cholmod.dll")

find_library(SUITESPARSE_LIB names suitesparseconfig
    HINTS ${CMAKE_LIBRARY_PATH} $ENV{SUITESPARSE_ROOT_LIB} "C:\\Users\\u0104126\\AppData\\Local\\Continuum\\miniconda3\\Library\\lib")

find_path(SUITESPARSE_INCLUDE names cholmod.h 
    HINTS ${CMAKE_INCLUDE_PATH} $ENV{SUITESPARSE_ROOT_INCLUDE})

message(STATUS "Root lib is set to: $ENV{SUITESPARSE_ROOT_LIB}")
message(STATUS "Root include is set to: $ENV{SUITESPARSE_ROOT_INCLUDE}")
message(STATUS ${CHOLMOD_LIB})
message(STATUS ${SUITESPARSE_INCLUDE})
message(STATUS ${SUITESPARSE_LIB})



# Add code coverage
option (COVERAGE "Perform code coverage" OFF)
if (COVERAGE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

# Add unit tests (if coverage is on, unit tests should also be on)
include(CMakeDependentOption)
CMAKE_DEPENDENT_OPTION(UNITTESTS "Perform unit testing" OFF
                       "NOT COVERAGE" ON)

message(STATUS "Unit testing is ${UNITTESTS}")
message(STATUS "Code coverage is ${COVERAGE}")

set(LIB_DIR ${CMAKE_CURRENT_LIST_DIR}/build/lib)
set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/build/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})

add_compile_options("-DDLONG")

#Rest of suitesparse
#add_subdirectory(cmakeSuitesparse)

add_library(qpalm SHARED
                src/newton.c
                src/util.c
                src/qpalm.c
                src/scaling.c
                src/linesearch.c
                src/lin_alg.c
                src/termination.c
                src/cholmod_interface.c
                src/nonconvex.c
                src/validate.c
                include/cholmod_interface.h
                include/global_opts.h
                include/constants.h
                include/nonconvex.h
                include/qpalm.h
                include/util.h
                include/scaling.h
                include/validate.h
                include/lin_alg.h
                include/linesearch.h
                include/types.h
                include/termination.h
                include/newton.h
                )

target_include_directories(qpalm
        PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                ${SUITESPARSE_INCLUDE})
        

target_link_libraries(qpalm ${CHOLMOD_LIB} ${SUITESPARSE_LIB})


# Add unittests
if (UNITTESTS)
    add_subdirectory(tests)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
