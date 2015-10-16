cmake_minimum_required(VERSION 3.1)
project(sssp)

set(CMAKE_C_STANDARD 99)
set(C_STANDARD_REQUIRED TRUE)

file(GLOB SRC_FILES *.c)

add_executable(sssp ${SRC_FILES})

install(TARGETS sssp DESTINATION bin)
