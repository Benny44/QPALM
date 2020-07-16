cmake_minimum_required(VERSION 3.0.2)
project(QPALM VERSION 0.1.0)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Windows-specific options
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

############################################################
################## Configure directories ###################
############################################################
set(LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/../lib)
set(INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../include)
make_directory(${INCLUDE_DIR})
if (NOT WIN32)
    execute_process(COMMAND bash -c "cp ${CMAKE_CURRENT_LIST_DIR}/include/* ${INCLUDE_DIR}")
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

############################################################
################### Compilation options ####################
############################################################
option (USE_CHOLMOD "Use CHOLMOD as the factorization package" OFF)
option (USE_LADEL "Use LADEL as the factorization package" OFF)
if (NOT USE_CHOLMOD AND NOT USE_LADEL)
    set (USE_LADEL ON)
    message (WARNING "No factorization package was selected. Using LADEL as the default.")
endif()

if (USE_CHOLMOD AND USE_LADEL)
    set (USE_CHOLMOD OFF)
    message (WARNING "Both CHOLMOD and LADEL selected. Using LADEL as the default.")
endif()

option (COVERAGE "Perform code coverage" OFF)
option (PYTHON "Compile for the python interface" OFF)
option (JULIA "Compile for the julia interface" OFF)
option (INTERFACES "Compile interfaces (qps/mtx)" ON)
option (PRINTING "Allow printing" ON)
option (PROF "Create executable demo for profiling with gprof" OFF)
option (PPROF "Create executable demo for profiling with cpuprofile" OFF)
option (TIMING "Enable timing" ON)
option (BUILD_SHARED_LIBS "Shared if on, static if off" ON)
# Add unit tests (if coverage is on, unit tests should also be on)
include(CMakeDependentOption)
CMAKE_DEPENDENT_OPTION(UNITTESTS "Perform unit testing" OFF
                       "NOT COVERAGE" ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDLONG")

if (USE_CHOLMOD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_CHOLMOD")
endif()

if (USE_LADEL)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_LADEL")
endif()

# Add code coverage
if (COVERAGE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

if (PYTHON)
    message("Compiling for the PYTHON interface")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPYTHON")
endif()

if (PRINTING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPRINTING")
endif()

if (PROF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
endif()

if (TIMING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPROFILING")
endif()

############################################################
###################### Find libraries ######################
############################################################
find_library(LAPACK names lapack
    HINTS $ENV{MINICONDA_LIB})
find_library(LAPACKE names lapacke
	HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})

if (USE_CHOLMOD)
    if (WIN32)
        find_library(CHOLMOD_LIB names cholmod 
            HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})

        find_library(SUITESPARSE_LIB names suitesparseconfig
            HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})

        find_path(SUITESPARSE_INCLUDE names cholmod.h 
            HINTS ${CMAKE_INCLUDE_PATH} $ENV{MINICONDA_INCLUDE}\\suitesparse)

        find_library(BLAS names openblas
            HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})
    else()
        find_library(MKL_INTEL_LP64 names mkl_intel_lp64 
            HINTS $ENV{MINICONDA_LIB})
        find_library(MKL_CORE names mkl_core 
            HINTS $ENV{MINICONDA_LIB})
        find_library(MKL_INTEL_THREAD names mkl_intel_thread 
            HINTS $ENV{MINICONDA_LIB})
        find_library(IOMP5 names iomp5 
            HINTS $ENV{MINICONDA_LIB})
        find_library(PTHREAD names pthread
            HINTS $ENV{MINICONDA_LIB})
    endif()
endif()

############################################################
################### Build the libraries ####################
############################################################
if (USE_CHOLMOD AND NOT WIN32)
    add_subdirectory(cmakeSuitesparse)
endif()

if (USE_LADEL)
    add_subdirectory(LADEL)
    if (NOT WIN32)
        execute_process(COMMAND bash -c "cp ${CMAKE_CURRENT_LIST_DIR}/LADEL/include/* ${INCLUDE_DIR}")
        execute_process(COMMAND bash -c "cp ${CMAKE_CURRENT_LIST_DIR}/LADEL/amd/Include/* ${INCLUDE_DIR}")
    endif()
endif()

if (UNIX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wconversion")
endif (UNIX)

SET(SRC 
        src/newton.c
        src/util.c
        src/qpalm.c
        src/scaling.c
        src/linesearch.c
        src/lin_alg.c
        src/termination.c
        src/solver_interface.c
        src/nonconvex.c
        src/validate.c
        src/iteration.c)

if (PYTHON)
    SET(SRC ${SRC} interfaces/python/qpalm_python.c)
endif()

if (JULIA)
    SET(SRC ${SRC} interfaces/QPALM.jl/qpalm_julia.c)
endif()

add_library(qpalm ${SRC})

############################################################
#################### Link the libraries ####################
############################################################
if (USE_CHOLMOD)
    if (NOT WIN32)
        target_include_directories(qpalm
            PUBLIC
                    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                    $ENV{MINICONDA_INCLUDE})
        target_link_libraries(qpalm ${MKL_INTEL_LP64} ${MKL_CORE} ${MKL_INTEL_THREAD} ${IOMP5} ${PTHREAD}
                            ${LAPACK} cholmod_new suitesparseconfig_new ${LAPACKE})
    else()
        target_include_directories(qpalm
            PUBLIC
                    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                    $ENV{MINICONDA_INCLUDE}
                    ${SUITESPARSE_INCLUDE})
        target_link_libraries(qpalm ${BLAS} ${LAPACK} ${CHOLMOD_LIB} ${SUITESPARSE_LIB} ${LAPACKE})
    endif()
endif()

if (USE_LADEL)
    target_include_directories(qpalm
            PUBLIC
                    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/LADEL/include>)
    target_link_libraries(qpalm ladel ${LAPACK} ${LAPACKE})
endif()


if (PYTHON)
    find_package(PythonInterp REQUIRED)
    # We don't really need the libs, but the headers.
    find_package(PythonLibs 3 REQUIRED)
    target_include_directories(qpalm PUBLIC ${PYTHON_INCLUDE_DIRS} interfaces/python/include)
    target_link_libraries(qpalm ${PYTHON_LIBRARIES})
endif()

############################################################
#################### Build executables #####################
############################################################
# Add unittests
if (UNITTESTS)
    add_subdirectory(tests)
    # Add testing
    include(CTest)
    enable_testing()
    add_test(NAME qpalm_tester COMMAND $<TARGET_FILE:run_all_tests>)
endif()

# Add interfaces
if (INTERFACES)
    option(PRINT_LATEX "Printout latex for qps" OFF)
    if (PRINT_LATEX)
        add_compile_options("-DPRINT_LATEX")
    endif()
    add_subdirectory(interfaces)
endif()

#Add profiling demo
if (PPROF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wl,--no-as-needed -lprofiler -Wl,--as-needed")
    add_subdirectory(profiling)
endif()

