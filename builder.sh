# !/bin/zsh

mkdir -p build
cd build
CC=gcc-12 CXX=g++-12 cmake -S .. -B . -DCMAKE_BUILD_TYPE=Debug
CC=gcc-12 CXX=g++-12 cmake --build . --parallel 8
