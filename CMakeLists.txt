cmake_minimum_required(VERSION 3.0.2)
project(QPALM VERSION 0.1.0)


set(LIB_DIR ${CMAKE_CURRENT_LIST_DIR}/build/lib)
set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/build/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})


# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)

make_directory(${INCLUDE_DIR})
message("Include dir is ${INCLUDE_DIR}")

if (NOT WIN32)
    execute_process(COMMAND bash -c "cp ${CMAKE_SOURCE_DIR}/include/* ${INCLUDE_DIR}")
endif()
#${CMAKE_SOURCE_DIR}/
# find_library(BLAS names openblas
#     HINTS /home/ben/miniconda3/lib)

message("Miniconda lib is $ENV{MINICONDA_LIB}")
message("Miniconda include is $ENV{MINICONDA_INCLUDE}")

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


find_library(LAPACK names lapack
    HINTS $ENV{MINICONDA_LIB})

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (WIN32)
    find_library(CHOLMOD_LIB names cholmod 
        HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})

    find_library(SUITESPARSE_LIB names suitesparseconfig
        HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})

    find_path(SUITESPARSE_INCLUDE names cholmod.h 
        HINTS ${CMAKE_INCLUDE_PATH} $ENV{MINICONDA_INCLUDE}\\suitesparse)

	find_library(BLAS names openblas
		HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})
endif()
# find_library(CHOLMOD_LIB names cholmod 
#     HINTS ${CMAKE_LIBRARY_PATH} $ENV{SUITESPARSE_ROOT_LIB})

# find_library(SUITESPARSE_LIB names suitesparseconfig
#     HINTS ${CMAKE_LIBRARY_PATH} $ENV{SUITESPARSE_ROOT_LIB})

# find_path(SUITESPARSE_INCLUDE names cholmod.h 
#     HINTS ${CMAKE_INCLUDE_PATH} $ENV{MINICONDA_INCLUDE})

# message("Suitesparse_include is ${SUITESPARSE_INCLUDE}")

find_library(LAPACKE names lapacke
	HINTS ${CMAKE_LIBRARY_PATH} $ENV{MINICONDA_LIB})


    #LAPACKE_INCLUDE assumed to be = to SUITESPARSE_INCLUDE
# find_path(LAPACK_INCLUDE names lapacke.h
#     HINTS ${CMAKE_INCLUDE_PATH} $ENV{SUITESPARSE_ROOT_INCLUDE})


# message(STATUS "Root lib is set to: $ENV{SUITESPARSE_ROOT_LIB}")
# message(STATUS "Root include is set to: $ENV{SUITESPARSE_ROOT_INCLUDE}")
# message(STATUS ${CHOLMOD_LIB})
# message(STATUS ${SUITESPARSE_INCLUDE})
# message(STATUS ${SUITESPARSE_LIB})


# Add code coverage
option (COVERAGE "Perform code coverage" OFF)
if (COVERAGE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

if (UNIX)
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wconversion")
endif (UNIX)

# Add unit tests (if coverage is on, unit tests should also be on)
include(CMakeDependentOption)
CMAKE_DEPENDENT_OPTION(UNITTESTS "Perform unit testing" OFF
                       "NOT COVERAGE" ON)

message(STATUS "Unit testing is ${UNITTESTS}")
message(STATUS "Code coverage is ${COVERAGE}")



add_compile_options("-DDLONG")
option (PRINTING "Allow printing" ON)
if (PRINTING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPRINTING")
endif()

option (PROF "Create executable demo for profiling with gprof" OFF)
if (PROF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
endif()


option (TIMING "Enable timing" ON)
if (TIMING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPROFILING")
endif()

option (BUILD_SHARED_LIBS "Shared if on, static if off" ON)

if (NOT WIN32)
    add_subdirectory(cmakeSuitesparse)
endif()

add_library(qpalm 
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
                src/iteration.c
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
                include/iteration.h
                )


                # ${SUITESPARSE_INCLUDE})
        
# find_package(OpenMP)
# if (OPENMP_FOUND)
#         target_link_libraries(qpalm ${OpenMP_C_LIBRARIES})
# endif()
# target_link_libraries(qpalm ${BLAS} ${CHOLMOD_LIB} ${SUITESPARSE_LIB} ${LAPACK_LIB})
# target_link_libraries(qpalm ${BLAS} ${LAPACK} cholmod_new ${LAPACK_LIB})
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
# target_link_libraries(qpalm blas lapack cholmod_new ${LAPACK_LIB})



option (PPROF "Create executable demo for profiling with cpuprofile" OFF)
if (PPROF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wl,--no-as-needed -lprofiler -Wl,--as-needed")
    # target_link_libraries(qpalm profiler)
    add_subdirectory(profiling)
endif()

# Add unittests
if (UNITTESTS)
    add_subdirectory(tests)
    # Add testing
    include(CTest)
    enable_testing()
    add_test(NAME qpalm_tester COMMAND $<TARGET_FILE:run_all_tests>)
endif()

# Add interfaces
option (INTERFACES "Compile interfaces" ON)
if (INTERFACES)
    add_subdirectory(interfaces)
    # add_subdirectory(mtx)
    # add_subdirectory(qps)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


#execute_process(COMMAND bash -c "ls ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "ls ${LIB_DIR}")
