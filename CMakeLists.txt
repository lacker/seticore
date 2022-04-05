cmake_minimum_required(VERSION 3.10)

project(seticore)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(HighFive)

add_executable(seticore main.cpp)
target_link_libraries(seticore HighFive)
