cmake_minimum_required(VERSION 4.0)

project(MeshSolve C CXX Fortran)

# list(PREPEND CMAKE_MODULE_PATH "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/clapack")

find_package(OpenBLAS CONFIG REQUIRED)
find_package(LAPACK REQUIRED)

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

target_include_directories(meshsolve 
    PRIVATE 
        ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
target_link_libraries(meshsolve
    PRIVATE OpenBLAS::OpenBLAS ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS}
    )
    