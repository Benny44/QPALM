cmake_minimum_required(VERSION 3.0.0)
project(SUITESPARSE VERSION 0.1.0)

#cmake for AMD, CAMD, COLAMD, CCOLAMD, CHOLMOD, SUITESPARSECONFIG
set(SUITESPARSE_CONFIG_PATH ${CMAKE_SOURCE_DIR}/suitesparse/SuiteSparse_config)
set(AMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/AMD) 
set(CAMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/CAMD) 
set(COLAMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/COLAMD) 
set(CCOLAMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/CCOLAMD) 
set(CHOLMOD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/CHOLMOD) 

#AMD
add_library(amd SHARED
${AMD_PATH}/Source/amd_post_tree.c
                ${AMD_PATH}/Source/amd_1.c
                ${AMD_PATH}/Source/amdbar.f
                ${AMD_PATH}/Source/amd_control.c
                ${AMD_PATH}/Source/amd_defaults.c
                ${AMD_PATH}/Source/amd_preprocess.c
                ${AMD_PATH}/Source/amd_info.c
                ${AMD_PATH}/Source/amd_global.c
                ${AMD_PATH}/Source/amd_postorder.c
                ${AMD_PATH}/Source/amd.f
                ${AMD_PATH}/Source/amd_order.c
                ${AMD_PATH}/Source/amd_aat.c
                ${AMD_PATH}/Source/amd_2.c
                ${AMD_PATH}/Source/amd_dump.c
                ${AMD_PATH}/Source/amd_valid.c
                ${AMD_PATH}/Include/amd_internal.h
                ${AMD_PATH}/Include/amd.h
                ${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h)

include_directories(amd 
        PUBLIC
                ${AMD_PATH}/Include
                ${SUITESPARSE_CONFIG_PATH})
                

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)