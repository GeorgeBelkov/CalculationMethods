cmake_minimum_required(VERSION 3.11)
project(cm_labs)
set(CMAKE_CXX_STANDARD 20)
include_directories(calc_met_lib lab1 lab2 lab3 lab4)
add_executable(lab lab3/main.cpp)