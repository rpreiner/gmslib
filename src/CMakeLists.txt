cmake_minimum_required(VERSION 3.15)

project(computeMixture)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(gmslib)

add_executable(computeMixture computeMixture.cpp)
target_link_libraries(computeMixture PRIVATE gmslib)
