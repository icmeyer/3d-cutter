cmake_minimum_required(VERSION 3.1...3.26)
project(
  3dcutter
  VERSION 0.1
  LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)

if(NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    message( FATAL_ERROR "3dcutter must be compiled with clang" )
endif()


# Library?
add_library(Lib3dcutter src/Cutting.cpp)

# Adding something we can run - Output name matches target name
add_executable(3dcutter src/Cutting.cpp)

target_link_libraries(3dcutter PRIVATE Lib3dcutter)
