#  cmakeLists.txt; 
#
#  Author(s): Ahmed Bachir <albachirou@gmail.com>
#             Kyoshe Winstone <wkyoshe@gmail.com>
#       Date: 2015-07-20
#
#  Université de Strasbourg.
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


#set the CMAKE_CXX_COMPILER
#set(CMAKE_CXX_COMPILER /usr/bin/clag++ 3.6.0)
cmake_minimum_required(VERSION 2.8.9 FATAL_ERROR)

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  cmake_minimum_required(VERSION 2.8)
  FIND_PATH(FEELPP_CMAKE_MODULES FindFeel++.cmake PATH  /usr/share/feel/cmake/modules/ /usr/local/share/feel/cmake/modules/ $ENV{FEELPP_DIR}/share/feel/cmake/modules )
  if ( FEELPP_CMAKE_MODULES )
    set(CMAKE_MODULE_PATH ${FEELPP_CMAKE_MODULES})
  else()
    message(FATAL_ERROR "Feel++ does not seem to have been installed on this platform")
  endif()
  Find_Package(Feel++)
  OPTION(FEELPP_ENABLE_SYSTEM_EIGEN3 "enable system eigen3 support" ON)
  if(ENV{CLN_DIR})
    INCLUDE_DIRECTORIES(${CLN_DIR}/include)
  endif()
endif()
SET(SRC
  02-environment
  03-mymesh
  04-myexpression
  05-myexporter
  07-myfunctionspace
  08-myintegrals
  #  model1.cpp
  #  model2.cpp
  #  myadvection.cpp
  #  myappwithoptions.cpp
  #  mybackend.cpp
  #  myfunctor.cpp
  #  mylaplacian.cpp
  #  mymeshfromstring.cpp
  #  mystokes.cpp
  )
foreach(app ${SRC})
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.cfg AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.geo )
    feelpp_add_application( ${app} SRCS ${app}.cpp  CFG ${app}.cfg GEO ${app}.geo)
  elseif ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.cfg )
    feelpp_add_application( ${app} SRCS ${app}.cpp  CFG ${app}.cfg )
  else()
    feelpp_add_application( ${app} SRCS ${app}.cpp )
  endif()
endforeach()
