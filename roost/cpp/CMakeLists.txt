cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(faust_anl)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

find_package(ROOT REQUIRED)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../target/cxxbridge/)

add_subdirectory("src")
