cmake_minimum_required(VERSION 3.0.0)
project(QPALM VERSION 0.1.0)

include(CTest)
enable_testing()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_path(SUITESPARSE_INCLUDE_DIR
    cholmod.h
HINTS ${CMAKE_CURRENT_LIST_DIR}/suitesparse/include)

find_library(SUITESPARSE_LIB
    names cholmod
    HINTS ${CMAKE_CURRENT_LIST_DIR}/suitesparse/lib)

add_library(qpalm SHARED src/newton.c
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
                ${SUITESPARSE_INCLUDE_DIR})
        
target_link_libraries(qpalm blas lapack ${SUITESPARSE_LIB})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
