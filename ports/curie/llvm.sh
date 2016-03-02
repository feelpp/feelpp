#! /bin/bash

cmake -G "Unix Makefiles" \
    -DCMAKE_INSTALL_PREFIX=/ccc/work/cont003/gen7335/chabanv/packages-install/llvm-3.7.0 \
    -DCMAKE_BUILD_TYPE=Release \
    -DGCC_INSTALL_PREFIX=/usr/local/gcc-4.9.1/ \
    -DCMAKE_C_COMPILER=/usr/local/gcc-4.9.1/bin/gcc \
    -DCMAKE_CXX_COMPILER=/usr/local/gcc-4.9.1/bin/g++ \
    -DBUILD_SHARED_LIBS=ON \
    -DLLVM_BUILD_TOOLS=OFF \
    -DCLANG_INCLUDE_DOCS=OFF -DCLANG_INCLUDE_TESTS=OFF \
    <path_llvm_src>
