cmake_minimum_required(VERSION 3.0.0)
project(QPALM_QPS VERSION 0.1.0)

add_executable(qpalm_qps 
                src/qpalm_qps.c
                src/qps_conversion.c
                src/index_hash.c)

message(${CMAKE_CURRENT_LIST_DIR}/include)
target_include_directories(qpalm_qps PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)

target_link_libraries(qpalm_qps qpalm)


