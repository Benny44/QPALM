cmake_minimum_required(VERSION 3.0.0)
project(SUITESPARSE VERSION 0.1.0)

# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)


# find_library(BLAS names mwblas
#     HINTS /home/ben/.MATLAB/R2015a/bin/glnxa64)

# find_library(LAPACK names mwlapack
#     HINTS /home/ben/.MATLAB/R2015a/bin/glnxa64)

# find_library(BLAS names openblas
#     HINTS /home/ben/miniconda3/lib)

# find_library(BLAS names mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 pthread
#     HINTS /home/ben/miniconda3/lib)

# find_library(LAPACK names lapack
#     HINTS /home/ben/miniconda3/lib)

#cmake for AMD, CAMD, COLAMD, CCOLAMD, CHOLMOD, SUITESPARSECONFIG
set(SUITESPARSE_CONFIG_PATH ${CMAKE_SOURCE_DIR}/suitesparse/SuiteSparse_config)
set(AMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/AMD) 
set(CAMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/CAMD) 
set(COLAMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/COLAMD) 
set(CCOLAMD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/CCOLAMD) 
set(CHOLMOD_PATH ${CMAKE_SOURCE_DIR}/suitesparse/CHOLMOD) 
set(METIS_PATH ${CMAKE_SOURCE_DIR}/suitesparse/metis-5.1.0)

#copy header files
# make_directory(${INCLUDE_DIR})
execute_process(COMMAND bash -c "cp ${AMD_PATH}/Include/amd.h ${INCLUDE_DIR}")
execute_process(COMMAND bash -c "cp ${CAMD_PATH}/Include/camd.h ${INCLUDE_DIR}")
execute_process(COMMAND bash -c "cp ${COLAMD_PATH}/Include/colamd.h ${INCLUDE_DIR}")
execute_process(COMMAND bash -c "cp ${CCOLAMD_PATH}/Include/ccolamd.h ${INCLUDE_DIR}")
execute_process(COMMAND bash -c "cp ${METIS_PATH}/include/metis.h ${INCLUDE_DIR}")
execute_process(COMMAND bash -c "cp ${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h ${INCLUDE_DIR}")
execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/* ${INCLUDE_DIR}")
execute_process(COMMAND bash -c "rm ${INCLUDE_DIR}/README.txt")

#execute_process(COMMAND bash -c "cp ${AMD_PATH}/Include/amd.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CAMD_PATH}/Include/camd.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${COLAMD_PATH}/Include/colamd.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CCOLAMD_PATH}/Include/ccolamd.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${METIS_PATH}/include/metis.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_blas.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_camd.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_check.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_cholesky.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_complexity.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_matrixops.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_config.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_modify.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_core.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_partition.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_function.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_supernodal.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_template.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_gpu.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_gpu_kernels.h ${INCLUDE_DIR}")
#execute_process(COMMAND bash -c "cp ${CHOLMOD_PATH}/Include/cholmod_io64.h ${INCLUDE_DIR}")  

add_compile_options("-DDLONG")

#SUITESPARSECONFIG
add_library(suitesparseconfig_new SHARED
${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.c
${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h
)

target_include_directories(suitesparseconfig_new
        PUBLIC 
                ${SUITESPARSE_CONFIG_PATH})

#AMD
add_library(amd_new SHARED
${AMD_PATH}/Source/amd_post_tree.c
                ${AMD_PATH}/Source/amd_1.c
                ${AMD_PATH}/Source/amd_control.c
                ${AMD_PATH}/Source/amd_defaults.c
                ${AMD_PATH}/Source/amd_preprocess.c
                ${AMD_PATH}/Source/amd_info.c
                ${AMD_PATH}/Source/amd_global.c
                ${AMD_PATH}/Source/amd_postorder.c
                ${AMD_PATH}/Source/amd_order.c
                ${AMD_PATH}/Source/amd_aat.c
                ${AMD_PATH}/Source/amd_2.c
                ${AMD_PATH}/Source/amd_dump.c
                ${AMD_PATH}/Source/amd_valid.c
                ${AMD_PATH}/Include/amd_internal.h
                ${AMD_PATH}/Include/amd.h
                ${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h)

target_include_directories(amd_new
        PUBLIC
                ${AMD_PATH}/Include
                ${SUITESPARSE_CONFIG_PATH})

target_link_libraries(amd_new suitesparseconfig_new)

#CAMD
add_library(camd_new SHARED
                ${CAMD_PATH}/Source/camd_dump.c
                ${CAMD_PATH}/Source/camd_order.c
                ${CAMD_PATH}/Source/camd_control.c
                ${CAMD_PATH}/Source/camd_2.c
                ${CAMD_PATH}/Source/camd_aat.c
                ${CAMD_PATH}/Source/camd_global.c
                ${CAMD_PATH}/Source/camd_info.c
                ${CAMD_PATH}/Source/camd_postorder.c
                ${CAMD_PATH}/Source/camd_1.c
                ${CAMD_PATH}/Source/camd_preprocess.c
                ${CAMD_PATH}/Source/camd_valid.c
                ${CAMD_PATH}/Source/camd_defaults.c
                ${CAMD_PATH}/Include/camd_internal.h
                ${CAMD_PATH}/Include/camd.h
                ${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h)

include_directories(camd_new
        PUBLIC
                ${CAMD_PATH}/Include
                ${SUITESPARSE_CONFIG_PATH})

target_link_libraries(camd_new suitesparseconfig_new)

#COLAMD
add_library(colamd_new SHARED
                ${COLAMD_PATH}/Source/colamd.c
                ${COLAMD_PATH}/Include/colamd.h
                ${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h
)

target_include_directories(colamd_new
        PUBLIC
                ${COLAMD_PATH}/Include
                ${SUITESPARSE_CONFIG_PATH})

target_link_libraries(colamd_new suitesparseconfig_new)

#CCOLAMD
add_library(ccolamd_new SHARED
                ${CCOLAMD_PATH}/Source/ccolamd.c
                ${CCOLAMD_PATH}/Include/ccolamd.h
                ${SUITESPARSE_CONFIG_PATH}/SuiteSparse_config.h
)

target_include_directories(ccolamd_new
        PUBLIC
                ${CCOLAMD_PATH}/Include
                ${SUITESPARSE_CONFIG_PATH})

target_link_libraries(ccolamd_new suitesparseconfig_new)

#CHOLMOD
add_library(cholmod_new SHARED
                ${CHOLMOD_PATH}/Core/cholmod_factor.c
                ${CHOLMOD_PATH}/Core/cholmod_band.c
                ${CHOLMOD_PATH}/Core/cholmod_version.c
                ${CHOLMOD_PATH}/Core/cholmod_transpose.c
                ${CHOLMOD_PATH}/Core/cholmod_common.c
                ${CHOLMOD_PATH}/Core/cholmod_dense.c
                ${CHOLMOD_PATH}/Core/cholmod_triplet.c
                ${CHOLMOD_PATH}/Core/cholmod_error.c
                ${CHOLMOD_PATH}/Core/cholmod_add.c
                ${CHOLMOD_PATH}/Core/cholmod_sparse.c
                ${CHOLMOD_PATH}/Core/cholmod_complex.c
                ${CHOLMOD_PATH}/Core/cholmod_change_factor.c
                ${CHOLMOD_PATH}/Core/cholmod_copy.c
                ${CHOLMOD_PATH}/Core/cholmod_memory.c
                ${CHOLMOD_PATH}/Core/cholmod_aat.c
                ${CHOLMOD_PATH}/Check/cholmod_check.c
                ${CHOLMOD_PATH}/Check/cholmod_write.c
                ${CHOLMOD_PATH}/Check/cholmod_read.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_analyze.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_rcond.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_solve.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_factorize.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_colamd.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_rowfac.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_etree.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_resymbol.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_amd.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_postorder.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_rowcolcounts.c
                ${CHOLMOD_PATH}/Cholesky/cholmod_spsolve.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_scale.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_submatrix.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_ssmult.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_drop.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_norm.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_symmetry.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_sdmult.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_vertcat.c
                ${CHOLMOD_PATH}/MatrixOps/cholmod_horzcat.c
                ${CHOLMOD_PATH}/Modify/cholmod_updown.c
                ${CHOLMOD_PATH}/Modify/cholmod_rowadd.c
                ${CHOLMOD_PATH}/Modify/cholmod_rowdel.c
                ${CHOLMOD_PATH}/Partition/cholmod_camd.c
                ${CHOLMOD_PATH}/Partition/cholmod_nesdis.c
                ${CHOLMOD_PATH}/Partition/cholmod_ccolamd.c
                ${CHOLMOD_PATH}/Partition/cholmod_csymamd.c
                ${CHOLMOD_PATH}/Partition/cholmod_metis.c
                ${CHOLMOD_PATH}/Supernodal/cholmod_super_numeric.c
                ${CHOLMOD_PATH}/Supernodal/cholmod_super_symbolic.c
                ${CHOLMOD_PATH}/Supernodal/cholmod_super_solve.c
                ${CHOLMOD_PATH}/Include/cholmod_core.h
                ${CHOLMOD_PATH}/Include/cholmod_modify.h
                ${CHOLMOD_PATH}/Include/cholmod_camd.h
                ${CHOLMOD_PATH}/Include/cholmod_gpu.h
                ${CHOLMOD_PATH}/Include/cholmod_gpu_kernels.h
                ${CHOLMOD_PATH}/Include/cholmod_function.h
                ${CHOLMOD_PATH}/Include/cholmod_partition.h
                ${CHOLMOD_PATH}/Include/cholmod_supernodal.h
                ${CHOLMOD_PATH}/Include/cholmod_internal.h
                ${CHOLMOD_PATH}/Include/cholmod_config.h
                ${CHOLMOD_PATH}/Include/cholmod_cholesky.h
                ${CHOLMOD_PATH}/Include/cholmod_blas.h
                ${CHOLMOD_PATH}/Include/cholmod_io64.h
                ${CHOLMOD_PATH}/Include/cholmod_complexity.h
                ${CHOLMOD_PATH}/Include/cholmod_check.h
                ${CHOLMOD_PATH}/Include/cholmod.h
                ${CHOLMOD_PATH}/Include/cholmod_template.h
                ${CHOLMOD_PATH}/Include/cholmod_matrixops.h)

target_include_directories(cholmod_new
        PUBLIC
                ${CHOLMOD_PATH}/Include
                ${AMD_PATH}/Include
                ${CAMD_PATH}/Include
                ${COLAMD_PATH}/Include
                ${CCOLAMD_PATH}/Include
                ${SUITESPARSE_CONFIG_PATH}
                ${METIS_PATH}/include)

message("LIB_DIR is ${LIB_DIR}")
find_library(METIS_LIB
    names metis
    HINTS ${LIB_DIR})


# find_package(OpenMP)
# if (OPENMP_FOUND)
#         target_link_libraries(cholmod_new ${OpenMP_C_LIBRARIES})
# endif()

# target_link_libraries(cholmod_new 
#                         ${BLAS} ${LAPACK} 
#                         ${METIS_LIB}
#                         suitesparseconfig_new
#                         amd_new
#                         camd_new
#                         colamd_new
#                         ccolamd_new)
target_link_libraries(cholmod_new 
                        ${MKL_INTEL_LP64} ${MKL_CORE} ${MKL_INTEL_THREAD} ${IOMP5} ${PTHREAD} 
                        ${LAPACK} 
                        ${METIS_LIB}
                        suitesparseconfig_new
                        amd_new
                        camd_new
                        colamd_new
                        ccolamd_new)

# target_link_libraries(cholmod_new 
#                         blas lapack
#                         ${METIS_LIB}
#                         suitesparseconfig_new
#                         amd_new
#                         camd_new
#                         colamd_new
#                         ccolamd_new)

