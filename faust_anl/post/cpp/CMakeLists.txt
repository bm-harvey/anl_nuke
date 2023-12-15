cmake_minimum_required(VERSION 3.16)
project("faust_anl")
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

find_package(ROOT REQUIRED)
set(CMAKE_CXX_STANDARD 14)



set(INC_DIR $ENV{HOME}/include)
set(LIB_DIR $ENV{HOME}/lib)

include_directories(${INC_DIR})

add_subdirectory("src")
