cmake_minimum_required(VERSION 3.10)

project(seticore)
set(CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

add_subdirectory(fmt)
add_subdirectory(HighFive)

add_executable(seticore main.cpp)
target_link_libraries(seticore fmt HighFive)
