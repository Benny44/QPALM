cmake_minimum_required(VERSION 3.0.0)
project(QPALMtest VERSION 0.1.0)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../build/test)

add_executable(run_all_tests
            run_all_tests.c
            src/test_dua_inf_qp.c
            src/test_validate.c
            src/test_update.c
            src/test_prim_inf_qp.c
            src/test_lin_alg.c
            src/test_basic_qp.c
            src/test_cholmod_interface.c
            src/test_nonconvex_qp.c
            src/test_degen_hess.c
            src/test_error_handling.c)

target_include_directories(run_all_tests 
        PRIVATE 
                include)

target_link_libraries(run_all_tests qpalm)
