cmake_minimum_required(VERSION 3.1...3.26)
project(
  cutter3d
  VERSION 0.1
  LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-L/System/Library/Frameworks -framework GLUT -framework OpenGL")

add_executable(cutter3d src/main.cpp)
# add_library(libcutter3d STATIC src/Cutting.h src/Cutting.cpp)
# target_link_libraries(cutter3d PRIVATE libcutter3d)

# INSTALL(TARGETS cutter3d
#     LIBRARY DESTINATION lib
#     ARCHIVE DESTINATION lib
#     )
# INSTALL(FILES src/Cutting.h DESTINATION include)
