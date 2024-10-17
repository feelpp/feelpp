# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#





SET( SRC
  marker.cpp
  entities.cpp
  regiontree.cpp kdtree.cpp
  functors.cpp #bvh.cpp
  meshbase.cpp
  remesher.cpp
  )


#if(TRUE)
#   find_package(rocthrust REQUIRED CONFIG PATHS "/nvme0/lemoinep/soft/rocThrust")
#   if (rocthrust_FOUND)
#      message("{rocthrust_LIBRARIES}") 
#      #set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/nvme0/lemoinep/soft/rocThrust/include")
#      set(CMAKE_PREFIX_PATH "/nvme0/lemoinep/soft/rocThrust/include")
#      list(APPEND include_dirs "/nvme0/lemoinep/soft/rocThrust/include")
#   else()
#      message(FATAL_ERROR "Could not find rocThrust")     
#   endif()
#endif()


#set(KERNEL_SOURCES
#   bvh.cpp
#)

#set_source_files_properties(${KERNEL_SOURCES}
#    PROPERTIES
#    HIP_SOURCE_PROPERTY_FORMAT 1
#    HIP_SEPARABLE_COMPILATION ON
#)



#find_package(HIP)
#enable_language(HIP)
#add_compile_definitions( __HIP_PLATFORM_AMD__)
#set_source_files_properties(bvh.hpp PROPERTIES LANGUAGE ${GPU_RUNTIME})
#set_target_properties(bvh.cpp PROPERTIES HIP_SEPARABLE_COMPILATION ON)
#set_source_files_properties(${KERNEL_SOURCES} PROPERTIES LANGUAGE ${GPU_RUNTIME})

#add_definitions(-DTHRUST_DEVICE_SYSTEM=${THRUST_DEVICE_SYSTEM})
#            include_directories(${rocTHRUST_INCLUDE_DIR})
#            link_directories(${rocTHRUST_INCLUDE_DIR})




file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp detail/*.hpp)

feel_append_src( feelmesh "${SRC};${HDR}")
