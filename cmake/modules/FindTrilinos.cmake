# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-12-21
#
#  Copyright (C) 2009 Université Joseph Fourier
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
FIND_PACKAGE(EXPAT)

FIND_PATH(TRILINOS_INCLUDE_DIR
    Teuchos_Utils.hpp
    PATHS /opt/trilinos/include/ /usr/include/trilinos /usr/ljk/include/trilinos
    $ENV{TRILINOS_DIR} 
    PATH_SUFFIXES
      include
    DOC "Directory where Trilinos header files are stored" )
MARK_AS_ADVANCED( TRILINOS_INCLUDE_DIR )

if( TRILINOS_INCLUDE_DIR )

  FIND_LIBRARY(TRILINOS_LIB_TEUCHOS   NAMES trilinos_teuchos teuchos     PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_EPETRA    NAMES trilinos_epetra epetra       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_EPETRAEXT NAMES trilinos_epetraext epetraext PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_TRIUTILS  NAMES trilinos_triutils triutils   PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_AZTECOO   NAMES trilinos_aztecoo aztecoo     PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_AMESOS    NAMES trilinos_amesos amesos       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_IFPACK    NAMES trilinos_ifpack ifpack       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_ML        NAMES trilinos_ml ml               PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_GALERI    NAMES trilinos_galeri galeri       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_NOX       NAMES trilinos_nox nox             PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
  FIND_LIBRARY(TRILINOS_LIB_NOXEPETRA NAMES trilinos_noxepetra noxepetra PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)
#  FIND_LIBRARY(TRILINOS_LIB_NOXLAPACK NAMES trilinos_noxlapack noxlapack PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib $ENV{TRILINOS_DIR}/lib)

  SET(TRILINOS_LIBRARIES
    ${TRILINOS_LIB_NOX}
    ${TRILINOS_LIB_NOXEPETRA}
#    ${TRILINOS_LIB_NOXLAPACK}
    ${TRILINOS_LIB_ML}
    ${TRILINOS_LIB_GALERI}
    ${TRILINOS_LIB_IFPACK}
    ${TRILINOS_LIB_AMESOS}
    ${TRILINOS_LIB_AZTECOO}
    ${TRILINOS_LIB_TRIUTILS}
    ${TRILINOS_LIB_EPETRAEXT}
    ${TRILINOS_LIB_EPETRA}
    ${TRILINOS_LIB_TEUCHOS}
    ${EXPAT_LIBRARIES}
    )

  SET(CMAKE_REQUIRED_FLAGS "${TRILINOS_LIBRARIES};${CMAKE_REQUIRED_FLAGS}")
  SET(CMAKE_REQUIRED_INCLUDES "${TRILINOS_INCLUDE_DIR};${CMAKE_REQUIRED_INCLUDES}")
  MESSAGE( STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")


  CHECK_INCLUDE_FILE_CXX(Teuchos_Utils.hpp FEELPP_HAS_TRILINOS_TEUCHOS )
  CHECK_INCLUDE_FILE_CXX(Epetra_Vector.h FEELPP_HAS_TRILINOS_EPETRA )
  CHECK_INCLUDE_FILE_CXX(EpetraExt_VectorIn.h FEELPP_HAS_TRILINOS_EPETRAEXT )
  CHECK_INCLUDE_FILE_CXX(Trilinos_Util_CrsMatrixGallery.h FEELPP_HAS_TRILINOS_TRIUTILS )
  CHECK_INCLUDE_FILE_CXX(AztecOO.h FEELPP_HAS_TRILINOS_AZTECOO )
  CHECK_INCLUDE_FILE_CXX(AztecOO.h FEELPP_HAS_AZTECOO_TEUCHOS )
  CHECK_INCLUDE_FILE_CXX(Amesos.h FEELPP_HAS_TRILINOS_AMESOS )
  CHECK_INCLUDE_FILE_CXX(Ifpack.h FEELPP_HAS_TRILINOS_IFPACK )
  CHECK_INCLUDE_FILE_CXX(ml_vec.h FEELPP_HAS_TRILINOS_ML )
  CHECK_INCLUDE_FILE_CXX(NOX_Utils.H FEELPP_HAS_TRILINOS_NOX )

  IF( FEELPP_HAS_TRILINOS_TEUCHOS )
    ADD_DEFINITIONS( -DFEELPP_HAS_TRILINOS -DFEELPP_HAS_TRILINOS_TEUCHOS -DFEELPP_HAS_TRILINOS_EPETRA -DFEELPP_HAS_TRILINOS_EPETRAEXT -DAVE_TRILINOS_TRIUTILS -DFEELPP_HAS_TRILINOS_AZTECOO -DFEELPP_HAS_AZTECOO_TEUCHOS -DFEELPP_HAS_TRILINOS_AMESOS -DFEELPP_HAS_TRILINOS_IFPACK  -DFEELPP_HAS_TRILINOS_ML -DML_MPI -DFEELPP_HAS_ML_TEUCHOS -DFEELPP_HAS_ML_EPETRA -DFEELPP_HAS_ML_AZTECOO )
  ENDIF( FEELPP_HAS_TRILINOS_TEUCHOS )


  MARK_AS_ADVANCED(
    TRILINOS_LIB_TEUCHOS
    TRILINOS_LIB_EPETRA
    TRILINOS_LIB_EPETRAEXT
    TRILINOS_LIB_TRIUTILS
    TRILINOS_LIB_AZTECOO
    TRILINOS_LIB_AMESOS
    TRILINOS_LIB_IFPACK
    TRILINOS_LIB_ML
    TRILINOS_LIB_GALERI
    TRILINOS_LIB_NOX
    TRILINOS_LIB_NOXEPETRA
#    TRILINOS_LIB_NOXLAPACK
    )


if ( TRILINOS_LIB_NOX AND
    TRILINOS_LIB_NOXEPETRA AND
    TRILINOS_LIB_ML AND
    TRILINOS_LIB_GALERI AND
    TRILINOS_LIB_IFPACK AND
    TRILINOS_LIB_AMESOS AND
    TRILINOS_LIB_AZTECOO AND
    TRILINOS_LIB_TRIUTILS AND
    TRILINOS_LIB_EPETRAEXT AND
    TRILINOS_LIB_EPETRA AND
    TRILINOS_LIB_TEUCHOS AND
    EXPAT_LIBRARIES )

  SET( TRILINOS_FOUND "YES")
endif ()

endif( TRILINOS_INCLUDE_DIR )
