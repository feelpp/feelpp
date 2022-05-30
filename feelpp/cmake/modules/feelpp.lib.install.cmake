###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-05-22
#
#  Copyright (C) 2014-2015 Feel++ Consortium
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
#
# Installing
#
INCLUDE(InstallRequiredSystemLibraries)

INSTALL(FILES ${CMAKE_SOURCE_DIR}/feelpp.version.cmake DESTINATION share/feelpp/feel/cmake/modules COMPONENT Devel)

feelpp_list_subdirs(feeldirs ${CMAKE_CURRENT_SOURCE_DIR})

foreach(includedir ${feeldirs})
  FILE(GLOB files "${includedir}/*.h*" )
  FILE(GLOB cppfiles "${includedir}/*.cpp" )
  INSTALL(FILES ${files} ${cppfiles} DESTINATION include/feelpp/feel/${includedir} COMPONENT Devel)
  
  feelpp_list_subdirs(feelsubdirs ${CMAKE_CURRENT_SOURCE_DIR}/${includedir})
  #message(STATUS "====== subdir ${includedir} : ${feelsubdirs}")
  foreach(subincludedir ${feelsubdirs})
    FILE(GLOB details "${includedir}/${subincludedir}/*.h*")
    #message(STATUS "subdir ${includedir}/${subincludedir}: ${details}")
    INSTALL(FILES ${details} DESTINATION include/feelpp/feel/${includedir}/${subincludedir}/ COMPONENT Devel)
  endforeach()
endforeach()
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/feel COMPONENT Devel)

# install cln headers
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/cln/include/cln/*.h")
INSTALL(FILES ${files} DESTINATION include/feelpp/cln COMPONENT Devel)
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/cln/lib/lib*" "${CMAKE_BINARY_DIR}/contrib/cln/lib64/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Libs)

# install gflags headers
# FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/gflags/include/gflags/*.h")
# INSTALL(FILES ${files} DESTINATION include/feel/gflags COMPONENT Devel)
# FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/gflags/lib/lib*" "${CMAKE_BINARY_DIR}/contrib/gflags/lib64/lib*")
# INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Libs)

# install glog headers
# FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/glog/include/glog/*.h")
# INSTALL(FILES ${files} DESTINATION include/feel/glog COMPONENT Devel)
# FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/glog/lib/lib*" "${CMAKE_BINARY_DIR}/contrib/glog/lib64/lib*")
# INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Libs)

# install cereal headers
if(FEELPP_HAS_CEREAL)
    FILE(GLOB_RECURSE files "${CMAKE_SOURCE_DIR}/contrib/cereal/include/*")
    FOREACH(fl IN LISTS files)
        string(REGEX REPLACE "${CMAKE_SOURCE_DIR}/contrib/cereal/include" "include/feel" fl1 ${fl})
        get_filename_component(dir ${fl1} DIRECTORY)
        INSTALL(FILES ${fl1} DESTINATION include/feelpp/cereal COMPONENT Devel)
    ENDFOREACH()
endif()

# feel++ config headers
FILE(GLOB files "${CMAKE_BINARY_DIR}/feelpp/feel/*.h")
INSTALL(FILES ${files} DESTINATION include/feelpp/feel COMPONENT Devel)

# feel++ precompiled headers.
if( FEELPP_ENABLE_PCH )
    file(GLOB files
        "${CMAKE_BINARY_DIR}/feelpp/feel/cotire/*.pch"
        "${CMAKE_BINARY_DIR}/feelpp/feel/cotire/*.cxx"
        "${CMAKE_BINARY_DIR}/feelpp/feel/cotire/*.hxx"
        )
    foreach(f IN LISTS files)
        install( FILES ${f} DESTINATION include/feelpp/feel/cotire COMPONENT Devel)
    endforeach()
endif()

