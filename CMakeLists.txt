cmake_minimum_required(VERSION 4.0)

project(MeshSolve C CXX Fortran)

find_package(OpenBLAS REQUIRED)

if(COMPILE_TO_DLL)
    add_library(meshsolve SHARED)
    set_target_properties(meshsolve PROPERTIES PREFIX "")

    target_sources(meshsolve 
        PRIVATE ${SRC_DIR}/meshsolve.c
        PUBLIC ${SRC_DIR}/meshsolve.h
    )
    message(STATUS "Compiling code to DLL library.")
else()
    add_executable(meshsolve ${SRC_DIR}/meshsolve.c)
    message(STATUS "Compiling code to executable.")
endif()

target_link_libraries(meshsolve
    PRIVATE OpenBLAS::OpenBLAS
    )
    