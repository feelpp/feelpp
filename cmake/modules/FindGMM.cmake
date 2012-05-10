# -*- mode: cmake -*-
#
#  This file is part of the Feel++ library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#             Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
#       Date: 2012-05-10
#
#  Copyright (C) 2012 Université Joseph Fourier
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
#  - GMM_INCLUDE_DIR
#  - GMM_FOUND

INCLUDE(CheckIncludeFileCXX)
# CHECK_INCLUDE_FILE_CXX(gmm/ANN.h FEELPP_HAS_GMM_H)


FIND_PATH(GMM_INCLUDE_DIR
  gmm.h
  PATHS /usr/include/ /usr/include/gmm /opt/local/include/gmm /usr/local/include/gmm  $ENV{GMM_DIR}/include/gmm
  DOC "Directory where gmm header files are stored" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMM "Could not find GMM headers" GMM_INCLUDE_DIR)
MARK_AS_ADVANCED(GMM_INCLUDE_DIR )
