###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 10 Jul 2017
#
#  Copyright (C) 2017 Feel++ Consortium
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
# MONGOCXX
#
option( FEELPP_ENABLE_MONGOCXX "Enable Mongocxx" ON )

if ( FEELPP_ENABLE_MONGOCXX )
  feelppContribPrepare( mongocxx )

  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${FEELPP_SOURCE_DIR}/contrib/mongocxx/cmake")
    find_package(LibBSON 1.5 )
    find_package(LibMongoC 1.5 )
    if ( LIBBSON_FOUND AND LIBMONGOC_FOUND )
      message(STATUS "[feelpp] found LibBSON and LibMongoC")
      set(BSONCXX_INLINE_NAMESPACE "")
      set(BSONCXX_HEADER_INSTALL_DIR "include/feelpp/bsoncxx/" CACHE INTERNAL "")

      #add_subdirectory(mongocxx)
      include_directories(${LIBBSON_INCLUDE_DIRS} ${LIBMONGOC_INCLUDE_DIRS})
      set(MONGOCXX_LIBRARIES "mongocxx_static ${LIBMONGOC_LIBRARIES} ${LIBBSON_LIBRARIES}")
      set(MONGOCXX_INCLUDE_DIRS "${LIBBSON_INCLUDE_DIRS} ${LIBMONGOC_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/src     ${CMAKE_CURRENT_BINARY_DIR}/src")
      set(LIBBSON_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
      list(APPEND LIBBSONCXX_INCLUDE_DIRS ${FEELPP_SOURCE_DIR}/contrib/mongocxx/src ${FEELPP_BINARY_DIR}/contrib/mongocxx/src)
      list(APPEND LIBMONGOCXX_INCLUDE_DIRS ${FEELPP_SOURCE_DIR}/contrib/mongocxx/src ${FEELPP_BINARY_DIR}/contrib/mongocxx/src)
      include_directories(${LIBMONGOCXX_INCLUDE_DIRS} ${LIBBSONCXX_INCLUDE_DIRS})
      set(FEELPP_HAS_BSONCXX 1)
      set(FEELPP_HAS_MONGOCXX 1)
      set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} MongoCXX" )
      #set(FEELPP_LIBRARIES "${MONGOCXX_LIBRARIES} ${FEELPP_LIBRARIES}")
      list(APPEND MONGO_LIBRARIES feelpp_mongocxx feelpp_bsoncxx ${LIBMONGOC_LIBRARIES} ${LIBBSON_LIBRARIES} )
      list(INSERT FEELPP_LIBRARIES 0 ${MONGO_LIBRARIES})
    endif()
  else()
      message(WARNING "MongoCXX was not found on your system. Either install it or set FEELPP_ENABLE_MONGOCXX to OFF.")
  endif()
endif()
