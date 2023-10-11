# !/bin/zsh

mkdir -p build
cd build
if [["$#" == 0]]; then
    CC=gcc-12 CXX=g++-12 cmake -S .. -B . -DCMAKE_BUILD_TYPE=Debug
elif [[ "$1" == "debug" ]]; then
    CC=gcc-12 CXX=g++-12 cmake -S .. -B . -DCMAKE_BUILD_TYPE=Debug
elif [[ "$1" == "release" ]]; then
    CC=gcc-12 CXX=g++-12 cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release
else
    echo "Invalid argument. Usage: ./builder.sh [debug|release]"
fi

CC=gcc-12 CXX=g++-12 cmake --build . --parallel 8
