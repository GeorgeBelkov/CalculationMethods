cmake_minimum_required(VERSION 3.11)
project(calc_met_lib)
set(CMAKE_CXX_STANDARD 20)
add_library(calc_met_lib STATIC Matrix.hpp Polinom.hpp)