cmake_minimum_required(VERSION 3.0)
project(quaternions)

message("System: " ${CMAKE_SYSTEM})
message("CMake:  " ${CMAKE_ROOT})
message("Using:  " ${CMAKE_CXX_COMPILER})

# boost - used for comparison in unit tests
find_package(Boost 1.5.9 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

set(SOURCE_FILES src/main.cpp src/Quaternion.h)

add_executable(quaternions ${SOURCE_FILES})
