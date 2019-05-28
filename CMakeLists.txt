cmake_minimum_required(VERSION 3.0.0)
project(QPALM VERSION 0.1.0)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(LIB_DIR ${CMAKE_CURRENT_LIST_DIR}/build/lib)
set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/build/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})

add_compile_options("-DDLONG")

#Rest of suitesparse
add_subdirectory(cmakeSuitesparse)

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

#message(${CMAKE_CURRENT_LIST_DIR}/suitesparse/CHOLMOD/Include)

target_include_directories(qpalm
        PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>)
        
target_link_libraries(qpalm blas lapack cholmod)

IF(CMAKE_BUILD_TYPE MATCHES "debug")
    add_subdirectory(tests)
ENDIF(CMAKE_BUILD_TYPE MATCHES "debug")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
