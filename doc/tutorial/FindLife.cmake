# - Find Life
# This module looks for Life (Library for the Finite Element Method) support
# it will define the following values
#  LIFE_INCLUDE_DIR = where life/lifecore/life.hpp can be found
#  LIFE_LIBRARY    = the library to link in
SET(LIFE_FOUND 0)
FIND_PATH(LIFE_INCLUDE_DIR lifeconfig.h  PATHS /usr/lib/life/include /opt/life/include /usr/ljk/include/life /usr/local  )

FIND_LIBRARY(LIFE_LIBRARY        life        PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib )
FIND_LIBRARY(LIFECORE_LIBRARY    lifecore    PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib )
FIND_LIBRARY(LIFEALG_LIBRARY     lifealg     PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib )
FIND_LIBRARY(LIFEMESH_LIBRARY    lifemesh    PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib )
FIND_LIBRARY(LIFEDISCR_LIBRARY   lifediscr   PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib )
FIND_LIBRARY(LIFEFILTERS_LIBRARY lifefilters PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib )

FIND_LIBRARY(BOOSTMPI_LIBRARY boost_mpi  PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib  )
FIND_LIBRARY(BOOSTSERIALIZATION_LIBRARY boost_serialization  PATHS /usr/lib/life/lib /opt/life/lib /usr/ljk/lib NO_DEFAULT_PATH NO_SYSTEM_ENVIRONMENT_PATH)


SET(LIFE_LIBRARIES
   ${LIFE_LIBRARY}
   ${LIFEFILTERS_LIBRARY} ${LIFEDISCR_LIBRARY} ${LIFEMESH_LIBRARY}  ${LIFEALG_LIBRARY} ${LIFECORE_LIBRARY} ${TRILINOS_LIBRARIES} ${PETSC_LIBRARIES}
   ${BOOSTMPI_LIBRARY} ${BOOSTSERIALIZATION_LIBRARY}
   )

IF(LIFE_INCLUDE_DIR AND LIFE_LIBRARIES)
   SET(LIFE_FOUND 1)
ENDIF(LIFE_INCLUDE_DIR AND LIFE_LIBRARIES)

MARK_AS_ADVANCED(
  LIFE_LIBRARY
  LIFECORE_LIBRARY
  LIFEALG_LIBRARY
  LIFEMESH_LIBRARY
  LIFEDISCR_LIBRARY
  LIFEFILTERS_LIBRARY

  BOOSTMPI_LIBRARY
  BOOSTSERIALIZATION_LIBRARY
  )
