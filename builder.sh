# !/bin/zsh

if [[ "$1" == "debug" ]]; then
    mkdir -p build
    cd build
    CC=gcc-12 CXX=g++-12 cmake -S .. -B . -DCMAKE_BUILD_TYPE=Debug
    CC=gcc-12 CXX=g++-12 cmake --build . --parallel 8
elif [[ "$1" == "release" ]]; then
    mkdir -p build
    cd build
    CC=gcc-12 CXX=g++-12 cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release
    CC=gcc-12 CXX=g++-12 cmake --build . --parallel 8
elif [[ "$1" == "clean" ]]; then
    rm -rf build
else
    echo "Invalid argument. Usage: ./builder.sh [debug|release|clean]"
fi
