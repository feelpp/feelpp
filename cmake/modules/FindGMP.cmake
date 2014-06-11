# -*- mode: cmake -*-
#
#  This file is part of the Feel++ library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-05-25
#
#  Copyright (C) 2014 Feel++ Consortium
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3.0 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# this files defines
#  - GMP_INCLUDE_DIR
#  - GMP_LIBRARIES
#  - GMP_FOUND

INCLUDE(CheckIncludeFileCXX)



FIND_LIBRARY( GMP_LIB gmp PATHS $ENV{GMP_DIR}/lib $ENV{GMP_LIB} NO_DEFAULT_PATH)
#FIND_LIBRARY( GMP_LIB gmp )
SET(GMP_LIBRARIES ${GMP_LIB} )

FIND_PATH(GMP_INCLUDE_DIR
  gmp.h
  PATHS  $ENV{GMP_DIR}/include $ENV{GMP_INC}
  DOC "Directory where GMP header files are stored"
  NO_DEFAULT_PATH  )

FIND_PATH(GMP_INCLUDE_DIR
  gmp.h
  DOC "Directory where GMP header files are stored" )

#CHECK_INCLUDE_FILE_CXX(gmp.h FEELPP_HAS_GMP_H)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP "Could not find GMP " GMP_INCLUDE_DIR GMP_LIBRARIES)
MARK_AS_ADVANCED(GMP_INCLUDE_DIR GMP_LIBRARIES )
