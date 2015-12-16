cmake_minimum_required(VERSION 3.0)
project(quaternions)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

message("System        : " ${CMAKE_SYSTEM})
message("CMake dir     : " ${CMAKE_ROOT})
message("Compiler      : " ${CMAKE_CXX_COMPILER})
message("Flags         : " ${CMAKE_CXX_FLAGS})

# boost - used for comparison in unit tests
find_package(Boost 1.5.9 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES src/main.cpp src/Quaternion.h)

add_executable(quaternions ${SOURCE_FILES})
