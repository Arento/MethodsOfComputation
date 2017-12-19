cmake_minimum_required(VERSION 3.9)
project(MethodsOfComputation)

set(CMAKE_CXX_STANDARD 17)

add_executable(Second Poly.cpp Poly.h InterpolationTools.h SecondMain.cpp)
add_executable(Third Poly.cpp Poly.h InterpolationTools.h ThirdMain.cpp)
add_executable(Fourth Poly.cpp Poly.h InterpolationTools.h FourthMain.cpp)
add_executable(Fifth Poly.cpp Poly.h InterpolationTools.h FifthMain.cpp)