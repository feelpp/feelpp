# \file Deps.cmake
# \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
# \date 2016-01-11

set(LLVM_VERSION "3.7")
set(LLVM_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/LLVM/${LLVM_VERSION})

set(OPENMPI_VERSION "1.10.1")
if(LOCAL_FILE_PATH)
    set(OPENMPI_URL "${LOCAL_FILE_PATH}/openmpi-${OPENMPI_VERSION}.tar.gz")
else()
    set(OPENMPI_URL "http://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-${OPENMPI_VERSION}.tar.gz")
endif()
set(OPENMPI_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/OpenMPI/${OPENMPI_VERSION})

set(GMSH_VERSION "2.11.0")
if(LOCAL_FILE_PATH)
    set(GMSH_URL "${LOCAL_FILE_PATH}/gmsh-${GMSH_VERSION}-source.tgz")
else()
    set(GMSH_URL "http://gmsh.info/src/gmsh-${GMSH_VERSION}-source.tgz")
endif()
set(GMSH_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/Gmsh/${GMSH_VERSION})

set(BOOST_VERSION "1.59.0")
string(REPLACE "." "_" BOOST_VERSION_US ${BOOST_VERSION})
if(LOCAL_FILE_PATH)
    set(BOOST_URL "${LOCAL_FILE_PATH}/boost_${BOOST_VERSION_US}.tar.gz")
else()
    set(BOOST_URL "http://sourceforge.net/projects/boost/files/boost/${BOOST_VERSION}/boost_${BOOST_VERSION_US}.tar.gz/download")
endif()
set(BOOST_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/Boost/${BOOST_VERSION})

set(PETSC_VERSION "3.6.3")
if(LOCAL_FILE_PATH)
    set(PETSC_URL "${LOCAL_FILE_PATH}/petsc-${PETSC_VERSION}.tar.gz")
else()
    set(PETSC_URL "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz")
endif()
set(PETSC_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/PETSc/${PETSC_VERSION})

set(SLEPC_VERSION "3.6.2")
if(LOCAL_FILE_PATH)
    set(SLEPC_URL "${LOCAL_FILE_PATH}/slepc-${SLEPC_VERSION}.tar.gz")
else()
    set(SLEPC_URL "http://slepc.upv.es/download/download.php?filename=slepc-${SLEPC_VERSION}.tar.gz")
endif()
set(SLEPC_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/SLEPc/${SLEPC_VERSION})

set(PARAVIEW_VERSION "4.4.0")
if(LOCAL_FILE_PATH)
    set(PARAVIEW_URL "${LOCAL_FILE_PATH}/ParaView-v4.4.0-source.tar.gz")
else()
    set(PARAVIEW_URL "http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.4&type=source&os=all&downloadFile=ParaView-v4.4.0-source.tar.gz")
endif()
set(PARAVIEW_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/ParaView/${PARAVIEW_VERSION})

set(_PATH "${OPENMPI_INSTALL_PREFIX}/bin:${GMSH_INSTALL_PREFIX}/bin:${PARAVIEW_INSTALL_PREFIX}/bin")
SET(_LD_LIBRARY_PATH "${OPENMPI_INSTALL_PREFIX}/lib:${GMSH_INSTALL_PREFIX}/lib:${BOOST_INSTALL_PREFIX}/lib:${PETSC_INSTALL_PREFIX}/lib:${PARAVIEW_INSTALL_PREFIX}/lib")

set(ENV{PATH} "${_PATH}:$ENV{PATH}")
set(ENV{LD_LIBRARY_PATH} "${_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
set(ENV{PETSC_DIR} "")
set(ENV{SLEPC_DIR} "")

message("$ENV{PATH}")
message("$ENV{LD_LIBRARY_PATH}")
message("$ENV{PETSC_DIR}")
message("$ENV{SLEPC_DIR}")

