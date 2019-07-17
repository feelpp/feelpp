#! /bin/bash

# generate the list of files to be added to the CMakeLists.txt
find geo -name "*.[meshdgo]*" | sort > list
