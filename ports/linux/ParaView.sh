#!/bin/bash

# Building ParaView 5.0.1

# For the display export over SSH to work, we need to use the old OpenGL backend
# See http://vtk.1045678.n5.nabble.com/GLXBadFBConfig-on-switching-from-VTK6-to-VTK7-td5738061.html
cmake $1 -DBUILD_TESTING=OFF -DVTK_RENDERING_BACKEND=OpenGL -DPARAVIEW_ENABLE_CATALYST=ON -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_USE_MPI=ON -DCMAKE_INSTALL_PREFIX=$2
make -j $3
