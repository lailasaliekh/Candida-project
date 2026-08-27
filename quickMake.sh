# Create a 
#

#!/bin/bash
NUM_THREADS=$(nproc --all)

echo "== Clearing old build files =="
rm -rf ./build

echo "== Creating fresh build directory =="
mkdir build

echo "== cd to build dir =="
cd ./build

echo "== Preparing MakeFile with CMake =="
cmake ..

echo "== Making =="
make -j $NUM_THREADS

echo "Run with the following command from the build directory: ./Main/main.out test/repeat0/ 0.6 1 1"
