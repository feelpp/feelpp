#! /bin/bash

brew tap homebrew/homebrew-science
brew update
brew install openmpi boost petsc slepc gmsh cln
mkdir build
cd build
../configure -r --cxxflags="-O2 -DNDEBUG"
make

