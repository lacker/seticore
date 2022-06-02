cmake_minimum_required(VERSION 3.10)
include(CMakePrintHelpers)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(seticore LANGUAGES CUDA C CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIR})
set(_hdf5_libs hdf5 hdf5_cpp)

find_package(Boost 1.40 COMPONENTS program_options REQUIRED)

add_subdirectory(fmt)

set(SOURCES
    dedoppler.cu
    dat_file_writer.cpp
    filterbank_file.cpp
    fil_file.cpp
    h5_file.cpp
    hit_file_writer.cpp
    recipe_file.cpp
    run_dedoppler.cpp
    )
               
add_executable(seticore ${SOURCES} main.cpp)
add_executable(test_beamforming ${SOURCES} test_beamforming.cpp)

set_property(TARGET seticore PROPERTY CUDA_STANDARD 14)
set_property(TARGET test_beamforming PROPERTY CUDA_STANDARD 14)

target_link_libraries(seticore fmt ${HDF5_C_LIBRARIES} Boost::program_options)
target_link_libraries(test_beamforming fmt ${HDF5_C_LIBRARIES} Boost::program_options)

