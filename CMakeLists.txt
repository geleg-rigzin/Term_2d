cmake_minimum_required(VERSION 3.20)

project(termo)

set(CMAKE_CXX_STANDART 11)
set(CMAKE_CXX_FLAGS "-O3")
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

set(MAIN_EXE main.cpp)

add_executable(main ${MAIN_EXE})
target_link_libraries(main LAPACK::LAPACK BLAS::BLAS)
