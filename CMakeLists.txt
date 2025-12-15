cmake_minimum_required(VERSION 4.0)

project(MeshSolve C CXX Fortran)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(OpenBLAS REQUIRED)

add_executable(meshsolve ${SRC_DIR}/meshsolve.c)

target_link_libraries(meshsolve
    PRIVATE OpenBLAS::OpenBLAS Python3::Python)
    