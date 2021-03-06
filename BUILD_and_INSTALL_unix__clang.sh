#!/bin/bash

echo

MYHOMEDEF=/mnt${HOME:1}/projects/share/cother

read -ep "Enter COTHER install path: " -i "${MYHOMEDEF}" MYHOME

echo
echo Install path: $MYHOME
echo

if [ ! -d build_clang ]; then mkdir build_clang || exit 1; fi
cd build_clang || exit 1


CC=clang CXX=clang++ \
cmake -DCMAKE_INSTALL_PREFIX=${MYHOME} \
    -DCMAKE_CUDA_FLAGS="-gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_75,code=compute_75" \
    ../src/  ||  (cd ..; exit 1)
#consider option: -DCMAKE_VERBOSE_MAKEFILE=ON

cmake --build . --config Release --target install  ||  (cd ..; exit 1)

cd ..

exit 0

