cmake_minimum_required(VERSION 3.10)
include(CMakePrintHelpers)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(seticore LANGUAGES CUDA C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIR})
set(_hdf5_libs hdf5 hdf5_cpp)

find_package(Boost 1.40 COMPONENTS program_options REQUIRED)

add_subdirectory(fmt)

add_executable(seticore main.cpp dedoppler.cu dat_file.cpp h5_file.cpp)
set_property(TARGET seticore PROPERTY CUDA_STANDARD 11)

target_link_libraries(seticore fmt ${HDF5_C_LIBRARIES} Boost::program_options)

