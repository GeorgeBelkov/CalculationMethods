cmake_minimum_required(VERSION 3.11)
project(lab1)
set(CMAKE_CXX_STANDARD 20)
target_include_directories(lab1 PUBLIC matrix)
target_link_libraries(lab1 matrix)
add_library(lab1 STATIC main.cpp)

