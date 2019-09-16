cmake_minimum_required(VERSION 3.0.0)
project(QPALMtest VERSION 0.1.0)

add_executable(demo_profile
            demo_profile.c
            load_data.c)

target_link_libraries(demo_profile qpalm)

if(UNIX)
target_link_libraries(demo_profile m)
endif(UNIX)

