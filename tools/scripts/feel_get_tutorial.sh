#!/bin/sh

mkdir -p $1
if [ -f /usr/share/doc/feel++-doc/examples/$1.cpp.gz ]; then zcat /usr/share/doc/feel++-doc/examples/$1.cpp > $1/$1.cpp; fi
if [ -f /usr/share/doc/feel++-doc/examples/$1.cpp ]; then cp /usr/share/doc/feel++-doc/examples/$1.cpp  $1/$1.cpp; fi
if [ -f /usr/share/doc/feel++-doc/examples/$1.hpp.gz ]; then zcat /usr/share/doc/feel++-doc/examples/$1.hpp > $1/$1.hpp; fi
if [ -f /usr/share/doc/feel++-doc/examples/$1.hpp ]; then cp /usr/share/doc/feel++-doc/examples/$1.hpp  $1/$1.hpp; fi
cat > $1/CMakeLists.txt << EOF
cmake_minimum_required(VERSION 2.8)

project($1 C CXX Fortran)

SET( CMAKE_MODULE_PATH  "/usr/share/feel/cmake/modules;\${CMAKE_CURRENT_SOURCE_DIR}" )

FIND_PACKAGE(Feel++ REQUIRED)

SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
include_directories(\${CMAKE_CURRENT_SOURCE_DIR})
add_executable($1 $1.cpp)
target_link_libraries($1 \${FEELPP_LIBRARIES})

message("")
message("You can now type 'make' to build the project \${PROJECT_NAME}")
message("")

EOF

echo "Now you can type the following commands to compile the tutorial $1"
echo " cd $1"
echo " cmake -DCMAKE_BUILD_TYPE=Release ."
echo " make"
