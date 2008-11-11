INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckLibraryExists)


SET( CMAKE_CXX_FLAGS "-pipe -Wall -O2 ")
SET( CMAKE_C_FLAGS "-pipe -Wall -O2")

project(life-tutorial)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})

#
# Blas and Lapack
#
FIND_LIBRARY(BLAS_LIBRARY  blas )
CHECK_LIBRARY_EXISTS("blas" "sdot_" "${BLAS_LIBRARY_DIR}" BLAS_LIBRARY_EXISTS)

FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack )
CHECK_LIBRARY_EXISTS("lapack" "dpotrf_" "${LAPACK_LIBRARY_DIR}" LAPACK_LIBRARY_EXISTS)
SET(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
LINK_DIRECTORIES(${LAPACK_LIBRARIES})

#
# Mpi
#
INCLUDE(FindMPI.cmake)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
LINK_DIRECTORIES(${MPI_LIBRARY_PATH})
ADD_DEFINITIONS( -DMPICH_IGNORE_CXX_SEEK )
SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")

#
# Vtk
#
FIND_PACKAGE(VTK QUIET)
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
LINK_DIRECTORIES(${VTK_LIBRARY_DIRS})

#
# Boost
#
FIND_PACKAGE(Boost)
#INCLUDE_DIRECTORIES(AFTER ${Boost_INCLUDE_DIRS})
SET(BOOST_LIBRARIES boost_program_options boost_filesystem )

#
# Metis
#
INCLUDE(FindMetis.cmake)
INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIR})
LINK_DIRECTORIES(${PARMETIS_LIBRARIES})

#
# Petsc
#
INCLUDE(FindPetsc.cmake)
INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR} ${PETSCCONF_INCLUDE_DIR})



#
# Trilinos
#
OPTION(ENABLE_TRILINOS "enable trilinos" OFF)
IF ( ENABLE_TRILINOS )
INCLUDE(FindTrilinos.cmake)
INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})
ENDIF( ENABLE_TRILINOS )

#
# Life
#
INCLUDE(FindLife.cmake)
INCLUDE_DIRECTORIES(${LIFE_INCLUDE_DIR})

SET(THELIBS
  ${LIFE_LIBRARIES}
  ${BOOST_LIBRARIES}
  ${PARMETIS_LIBRARIES}
  ${MPI_LIBRARIES}
  ${LAPACK_LIBRARIES}
  )

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})




add_executable (laplacian laplacian.cpp)
target_link_libraries(laplacian ${THELIBS} )

add_executable (myintegrals myintegrals.cpp)
target_link_libraries(myintegrals ${THELIBS} )

add_executable (mymesh mymesh.cpp)
target_link_libraries(mymesh ${THELIBS} )

#add_executable (mympiapp mympiapp.cpp)
#target_link_libraries(mympiapp  ${THELIBS})

add_executable (myapp myapp.cpp)
target_link_libraries(myapp  ${THELIBS})


#
# Latex
#
OPTION(ENABLE_LATEX "enable latex" OFF)
IF( ENABLE_LATEX )
INCLUDE(UseLATEX.cmake)
SET(LATEX_OUTPUT_PATH latex)
ADD_LATEX_DOCUMENT(
  tutorial.tex
  INPUTS myapp.cpp mympiapp.cpp mymesh.cpp myintegrals.cpp laplacian.cpp
  IMAGE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  IMAGES mymeshpartition.png
  BIBFILES tutorial.bib
  DEFAULT_PDF
  )
ENDIF( ENABLE_LATEX )

#
# testing
#
ENABLE_TESTING()
ADD_TEST(myapp myapp )
ADD_TEST(mymesh mymesh )
ADD_TEST(myintegrals myintegrals )
ADD_TEST(laplacian laplacian )

#
# Install
#
INSTALL(FILES README COPYING DESTINATION doc)
INSTALL(FILES ${CMAKE_BINARY_DIR}/myapp ${CMAKE_BINARY_DIR}/mympiapp ${CMAKE_BINARY_DIR}/mymesh ${CMAKE_BINARY_DIR}/myintegrals ${CMAKE_BINARY_DIR}/laplacian DESTINATION bin)


#
# Cpack
#
INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_GENERATOR "TGZ")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "life-tutorial")
SET(CPACK_PACKAGE_VENDOR "Christophe Prud'homme")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "9")
SET(CPACK_PACKAGE_VERSION_PATCH "8")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "life-tutorial")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "CPackSourceConfig.cmake")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "life-tutorial-0.9.8")
SET(CPACK_SOURCE_STRIP_FILES "")
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
set(CPACK_SOURCE_IGNORE_FILES
"~$"
"\\\\.cvsignore$"
"^${PROJECT_SOURCE_DIR}.*/.svn/"
"^${PROJECT_SOURCE_DIR}/CMakeFiles/"
"^${PROJECT_SOURCE_DIR}/latex/"
"^${PROJECT_SOURCE_DIR}/*~"
)


INCLUDE(CPack)
