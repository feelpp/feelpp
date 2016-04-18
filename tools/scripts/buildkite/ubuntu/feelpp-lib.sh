#!/bin/bash

set -e

echo '--- build dirctory'
mkdir build
cd build

echo '--- configure -r'
../configure -r

echo '--- make feelpp library'
make -j8 feelpp
