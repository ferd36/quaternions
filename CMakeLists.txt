cmake_minimum_required(VERSION 3.3)
project(quaternions)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES
    src/main.cpp)

add_executable(quaternions ${SOURCE_FILES})