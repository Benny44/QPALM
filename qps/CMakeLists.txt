cmake_minimum_required(VERSION 3.0.0)
project(QPALMqps VERSION 0.1.0)

add_executable(qpalm_qps qpalm_qps.c)

target_link_libraries(qpalm_qps qpalm)


