cmake_minimum_required(VERSION 3.11)
project(lab2)
set(CMAKE_CXX_STANDARD 20)
target_include_directories(lab2 PUBLIC matrix)
target_link_libraries(lab2 matrix)
add_library(lab2 STATIC main.cpp)