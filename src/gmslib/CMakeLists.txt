cmake_minimum_required(VERSION 3.15)

project(gmslib)

find_package(OpenMP REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} 
    INTERFACE 
        .
        ../
)
target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_C)
