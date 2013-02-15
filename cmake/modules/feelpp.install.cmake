###  coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2013-02-08
#
#  Copyright (C) 2013 Université de Strasbourg
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
set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEELPP_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEELPP_DATADIR ${CMAKE_INSTALL_PREFIX}/share/feel )
CONFIGURE_FILE(feelconfig.h.in feel/feelconfig.h  @ONLY)
CONFIGURE_FILE(feelinfo.h.in feel/feelinfo.h  @ONLY)

#
# Packaging
#
INCLUDE(InstallRequiredSystemLibraries)


foreach(includedir feelcore feelalg feelmesh feelpoly feelfilters feeldiscr feelvf feelmaterial feelsystem )
  FILE(GLOB files "feel/${includedir}/*.hpp")
  FILE(GLOB cppfiles "feel/${includedir}/*.cpp")
  INSTALL(FILES ${files} ${cppfiles} DESTINATION include/feel/${includedir} COMPONENT Devel)
endforeach()
FILE(GLOB files "feel/*.hpp")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)

# install cln headers
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/cln/include/cln/*.h")
INSTALL(FILES ${files} DESTINATION include/feel/cln COMPONENT Devel)
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/cln/lib/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Devel)

# install gflags headers
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/gflags/include/gflags/*.h")
INSTALL(FILES ${files} DESTINATION include/feel/gflags COMPONENT Devel)
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/gflags/lib/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Devel)

# install glog headers
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/glog/include/glog/*.h")
INSTALL(FILES ${files} DESTINATION include/feel/glog COMPONENT Devel)
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/glog/lib/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Devel)

# install matheval headers
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/libmatheval/include/*.h")
INSTALL(FILES ${files} DESTINATION include/feel/matheval COMPONENT Devel)
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/contrib/libmatheval/lib/lib*")
INSTALL(FILES ${files} DESTINATION lib/ COMPONENT Devel)

# gmm
IF ( NOT GMM_FOUND )
  FILE(GLOB files "contrib/gmm/include/*.h")
  INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)
ENDIF()
# feel++ config headers
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/feel/*.h")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)

#message(STATUS "apps: ${FEELPP_INSTALL_APPS}")
#INSTALL(FILES ${FEELPP_INSTALL_APPS} DESTINATION bin/ COMPONENT Bin)

IF(FEELPP_ENABLE_DOCUMENTATION)
# documentation and examples
#  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/feel_get_tutorial.sh DESTINATION bin COMPONENT Doc)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/manual/feel-manual.pdf DESTINATION share/doc/feel COMPONENT Doc)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/api/html DESTINATION share/doc/feel COMPONENT Doc
        PATTERN ".svn" EXCLUDE PATTERN ".git" EXCLUDE )
FILE(GLOB examples "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/*.*pp")
INSTALL(FILES ${examples} DESTINATION share/doc/feel/examples/ COMPONENT Doc)
ENDIF()
