# - Try to find Berkeley DB
# Once done this will define
#
#  BERKELEY_DB_FOUND - system has Berkeley DB
#  BERKELEY_DB_INCLUDE_DIR - the Berkeley DB include directory
#  BERKELEY_DB_LIBRARIES - Link these to use Berkeley DB
#  BERKELEY_DB_DEFINITIONS - Compiler switches required for using Berkeley DB

# Copyright (c) 2006, Alexander Dymo, <adymo@kdevelop.org>
# Copyright (c) 2011, Université Joseph Fourier
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

FIND_PATH(BERKELEY_DB_INCLUDE_DIR db_cxx.h
  /usr/include/db4
  /usr/local/include/db4
  /usr/include/db5
  /usr/local/include/db5
)

FIND_LIBRARY(BERKELEY_DB_LIBRARIES NAMES db_cxx )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Berkeley "Could not find Berkeley DB >= 4.1" BERKELEY_DB_INCLUDE_DIR BERKELEY_DB_LIBRARIES)
# show the BERKELEY_DB_INCLUDE_DIR and BERKELEY_DB_LIBRARIES variables only in the advanced view
MARK_AS_ADVANCED(BERKELEY_DB_INCLUDE_DIR BERKELEY_DB_LIBRARIES )
