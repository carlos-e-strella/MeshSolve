cmake_minimum_required(VERSION 4.0)

project(MeshSolve C CXX Fortran)

find_package(OpenBLAS REQUIRED)

add_executable(meshsolve ${SRC_DIR}/main.c)

target_link_libraries(meshsolve
    PRIVATE OpenBLAS::OpenBLAS)
    