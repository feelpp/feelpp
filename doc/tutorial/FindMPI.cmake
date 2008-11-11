# - Find MPI
# This module looks for MPI (Message Passing Interface) support
# it will define the following values
#  MPI_INCLUDE_PATH = where mpi.h can be found
#  MPI_LIBRARY    = the library to link in (mpi mpich etc)

FIND_PATH(MPI_INCLUDE_PATH mpi.h
          /usr/local/mpich/include
          /usr/local/include
          /usr/include
          /usr/include/mpi
          /usr/local/mpi/include
          "$ENV{ProgramFiles}/MPICH/SDK/Include"
          "$ENV{ProgramFiles}/MPICH2/include"
          "C:/Program Files/MPICH/SDK/Include"
)
FIND_PATH(MPI_LIBRARY_PATH
   NAMES libmpi.a libmpich.a libmpi.so libmpich.so
   PATHS
          /usr/local/mpich/lib
          /usr/local/lib
          /usr/lib
)
FIND_LIBRARY(PMPI_LIBRARY
             NAMES pmpich
             PATHS /usr/local/mpich/lib /usr/lib /usr/local/lib /usr/local/mpi/lib
             "$ENV{ProgramFiles}/MPICH/SDK/Lib"
             "$ENV{ProgramFiles}/MPICH2/Lib"
             "C:/Program Files/MPICH/SDK/Lib"
)
IF ( PMPI_LIBRARY )
  SET(MPI_EXTRA_LIBRARY_BEFORE ${PMPI_LIBRARY})
ENDIF( PMPI_LIBRARY )

FIND_LIBRARY(MPI_LIBRARY
             NAMES mpich2 mpi mpich
             PATHS /usr/local/mpich/lib /usr/lib /usr/local/lib /usr/local/mpi/lib
             "$ENV{ProgramFiles}/MPICH/SDK/Lib"
             "$ENV{ProgramFiles}/MPICH2/Lib"
             "C:/Program Files/MPICH/SDK/Lib"
)

FIND_LIBRARY(MPI_EXTRA_LIBRARY
             NAMES mpi++ mpicxx mpichcxx
             PATHS /usr/local/mpich/lib /usr/lib /usr/local/lib /usr/local/mpi/lib
             "$ENV{ProgramFiles}/MPICH/SDK/Lib"
             "C:/Program Files/MPICH/SDK/Lib"
             DOC "If a second mpi library is necessary, specify it here.")

SET(MPI_LIBRARIES ${MPI_EXTRA_LIBRARY} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY_BEFORE})
MARK_AS_ADVANCED(MPI_INCLUDE_PATH MPI_EXTRA_LIBRARY_BEFORE MPI_LIBRARY MPI_EXTRA_LIBRARY)
