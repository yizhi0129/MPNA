cmake_minimum_required(VERSION 3.12)
project(mpna-matrix-matrix C)

set(CMAKE_C_STANDARD 11)

add_executable(matrix-matrix main.c)
target_compile_options(matrix-matrix PRIVATE "-Wno-incompatible-pointer-types")
