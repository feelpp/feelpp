###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-05-22
#
#  Copyright (C) 2014 Feel++ Consortium
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
feelpp_list_subdirs(feeldirs ${CMAKE_CURRENT_SOURCE_DIR})

foreach(includedir ${feeldirs})
  FILE(GLOB files "${includedir}/*.hpp" )
  FILE(GLOB cppfiles "${includedir}/*.cpp" )
  INSTALL(FILES ${files} ${cppfiles} DESTINATION include/feel/${includedir} COMPONENT Devel)
  FILE(GLOB details "${includedir}/detail/*.hpp")
  INSTALL(FILES ${details} DESTINATION include/feel/${includedir}/detail COMPONENT Devel)
endforeach()
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)

# install cln headers
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/cln/include/cln/*.h")
INSTALL(FILES ${files} DESTINATION include/feel/cln COMPONENT Devel)
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/cln/lib/lib*" "${CMAKE_BINARY_DIR}/contrib/cln/lib64/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Devel)

# install gflags headers
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/gflags/include/gflags/*.h")
INSTALL(FILES ${files} DESTINATION include/feel/gflags COMPONENT Devel)
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/gflags/lib/lib*" "${CMAKE_BINARY_DIR}/contrib/gflags/lib64/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Devel)

# install glog headers
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/glog/include/glog/*.h")
INSTALL(FILES ${files} DESTINATION include/feel/glog COMPONENT Devel)
FILE(GLOB files "${CMAKE_BINARY_DIR}/contrib/glog/lib/lib*" "${CMAKE_BINARY_DIR}/contrib/glog/lib64/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Devel)

# install cereal headers
if(FEELPP_HAS_CEREAL)
    FILE(GLOB_RECURSE files "${CMAKE_SOURCE_DIR}/contrib/cereal/include/*")
    FOREACH(fl IN LISTS files)
        string(REGEX REPLACE "${CMAKE_SOURCE_DIR}/contrib/cereal/include" "include/feel" fl1 ${fl})
        get_filename_component(dir ${fl1} DIRECTORY)
        INSTALL(FILES ${fl1} DESTINATION include/feel/cereal COMPONENT Devel)
    ENDFOREACH()
endif()

# feel++ config headers
FILE(GLOB files "${CMAKE_BINARY_DIR}/feel/*.h")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)
