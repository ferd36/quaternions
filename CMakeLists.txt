cmake_minimum_required(VERSION 3.3)
project(quaternions)

message("System: " ${CMAKE_SYSTEM})
message("CMake:  " ${CMAKE_ROOT})
message("Using:  " ${CMAKE_CXX_COMPILER})

# debug
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG -g -O0 -Wall -Wextra -msse4.1 -DHAVE_CXX0X -march=native")

# profile
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse4.1 -DHAVE_CXX0X -march=native -gdwarf-3 -g -fprofile-arcs -ftest-coverage")

# release
# -Wcast-align
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -DNDEBUG -msse4.1 -DHAVE_CXX0X -march=native")

# boost
find_package(Boost 1.5.9 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse4.1 -march=native -std=c++14")

add_subdirectory(vectorclass)

set(SOURCE_FILES
    src/main.cpp src/Quaternion.h src/generators.h)

add_executable(quaternions ${SOURCE_FILES})
target_link_libraries(vectorclass)