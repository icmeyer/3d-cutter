cmake_minimum_required(VERSION 3.1...3.26)
project(
  cutter3d
  VERSION 0.1
  LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)

add_library(cutter3d STATIC src/Cutting.cpp)

add_executable(cutter3d_exec src/Cutting.cpp)

INSTALL(TARGETS cutter3d
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
INSTALL(FILES src/Cutting.h DESTINATION include)
