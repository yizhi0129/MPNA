cmake_minimum_required(VERSION 3.12)
project(mpna C CXX)


set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)

add_subdirectory(01-matrix-matrix)
add_subdirectory(02-iterative-relaxation)
