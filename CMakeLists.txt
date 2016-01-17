cmake_minimum_required(VERSION 3.0) # travis doesn't want 3.4 I'm using locally
project(quaternions)

set(CMAKE_CXX_FLAGS "-Wall -std=c++11")

message("System        : " ${CMAKE_SYSTEM})
message("CMake dir     : " ${CMAKE_ROOT})
message("Compiler      : " ${CMAKE_CXX_COMPILER})
message("Flags         : " ${CMAKE_CXX_FLAGS})

# boost - used for comparison in unit tests
# uncomment only if needing to use a specific location.
# in travis, it assumes boost is in /usr/include/boost
#find_package(Boost 1.5.9 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES test/unit_tests.cpp include/quaternion.h include/quaternion_utils.h include/quaternion_io.h)

add_executable(quaternions ${SOURCE_FILES})
