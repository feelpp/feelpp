# -*- mode: makefile -*-
#
#  This file is part of the Life library
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
find_package(EXPAT)

FIND_PATH(TRILINOS_INCLUDE_DIR Teuchos_Utils.hpp PATHS /opt/trilinos/include/ /usr/include/trilinos /usr/ljk/include/trilinos )
MARK_AS_ADVANCED( TRILINOS_INCLUDE_DIR )
if( TRILINOS_INCLUDE_DIR )

  SET(CMAKE_REQUIRED_INCLUDES "${TRILINOS_INCLUDE_DIR};${CMAKE_REQUIRED_INCLUDES}")
  MESSAGE( STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")


  CHECK_INCLUDE_FILE_CXX("Teuchos_Utils.hpp" HAVE_TRILINOS_TEUCHOS  )
  IF ( TRILINOS_LIB_TEUCHOS )
    SET( CMAKE_REQUIRED_FLAGS "-lteuchos;${CMAKE_REQUIRED_FLAGS}" )
  ENDIF( TRILINOS_LIB_TEUCHOS )

  CHECK_INCLUDE_FILE_CXX(Epetra_Vector.h HAVE_TRILINOS_EPETRA )
  CHECK_INCLUDE_FILE_CXX(EpetraExt_VectorIn.h HAVE_TRILINOS_EPETRAEXT )
  CHECK_INCLUDE_FILE_CXX(Trilinos_Util_CrsMatrixGallery.h HAVE_TRILINOS_TRIUTILS )
  CHECK_INCLUDE_FILE_CXX(AztecOO.h HAVE_TRILINOS_AZTECOO )
  CHECK_INCLUDE_FILE_CXX(AztecOO.h HAVE_AZTECOO_TEUCHOS )
  CHECK_INCLUDE_FILE_CXX(Amesos.h HAVE_TRILINOS_AMESOS "-ltrilinos_teuchos ${MPI_EXTRA_LIBRARY} ${MPI_LIBRARY}")
  CHECK_INCLUDE_FILE_CXX(Ifpack.h HAVE_TRILINOS_IFPACK "-ltrilinos_teuchos ${MPI_EXTRA_LIBRARY} ${MPI_LIBRARY}")
  CHECK_INCLUDE_FILE_CXX(ml_vec.h HAVE_TRILINOS_ML "-ltrilinos_teuchos ${MPI_EXTRA_LIBRARY} ${MPI_LIBRARY}")
  CHECK_INCLUDE_FILE_CXX(NOX_Utils.H HAVE_TRILINOS_NOX "-ltrilinos_nox -ltrilinos_noxepetra  -lteuchos ${MPI_EXTRA_LIBRARY} ${MPI_LIBRARY}")

  IF( HAVE_TRILINOS_TEUCHOS )
    ADD_DEFINITIONS( -DHAVE_TRILINOS -DHAVE_TRILINOS_TEUCHOS -DHAVE_TRILINOS_EPETRA -DHAVE_TRILINOS_EPETRAEXT -DAVE_TRILINOS_TRIUTILS -DHAVE_TRILINOS_AZTECOO -DHAVE_AZTECOO_TEUCHOS -DHAVE_TRILINOS_AMESOS -DHAVE_TRILINOS_IFPACK  -DHAVE_TRILINOS_ML -DML_MPI -DHAVE_ML_TEUCHOS -DHAVE_ML_EPETRA -DHAVE_ML_AZTECOO )
  ENDIF( HAVE_TRILINOS_TEUCHOS )

  FIND_LIBRARY(TRILINOS_LIB_TEUCHOS   NAMES trilinos_teuchos teuchos     PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_EPETRA    NAMES trilinos_epetra epetra       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_EPETRAEXT NAMES trilinos_epetraext epetraext PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_TRIUTILS  NAMES trilinos_triutils triutils   PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_AZTECOO   NAMES trilinos_aztecoo aztecoo     PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_AMESOS    NAMES trilinos_amesos amesos       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_IFPACK    NAMES trilinos_ifpack ifpack       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_ML        NAMES trilinos_ml ml               PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_GALERI    NAMES trilinos_galeri galeri       PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_NOX       NAMES trilinos_nox nox             PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
  FIND_LIBRARY(TRILINOS_LIB_NOXEPETRA NAMES trilinos_noxepetra noxepetra PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )
#  FIND_LIBRARY(TRILINOS_LIB_NOXLAPACK NAMES trilinos_noxlapack noxlapack PATHS /usr/lib /opt/trilinos/lib /usr/ljk/lib )

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

#SET( TRILINOS_FOUND )
endif()
