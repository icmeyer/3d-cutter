#!/bin/bash
# clang++ -o cutter main.cpp -L/System/Library/Frameworks -framework GLUT -framework OpenGL -std=c++11
# clang++ -o cutter Cutting.cpp -std=c++11
# clang++ -o cutter test.cpp -std=c++11
rm -rf build
mkdir build 
cd build
CXX=/usr/bin/clang++ cmake ../
make -j4
