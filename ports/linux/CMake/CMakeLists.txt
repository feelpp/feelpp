# \file CMakeLists.txt
# \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
# \date 2016-01-11

# This script allows to build Feel++ main dependencies with CMake
# It allows to build those dependencies in parallel with CMake dependency system
# To do so, first create a build directory and go in it, then type the following commands:
# - cmake <FEELPP_SRC_PATH>/ports/linux/CMake -DCMAKE_INSTALL_PREFIX=<path_where_to_install_deps>
# - make -j <NPROC>
# and you dependecies should be installed in subdirectories of <path_where_to_install_deps>

project(feelpp-dependencies)
cmake_minimum_required (VERSION 2.8.7)

include(ExternalProject)

include(Deps.cmake)

OPTION(ENABLE_BUILD_CLANG "Enable Clang compilation" OFF)
OPTION(ENABLE_GUI "Enable GUI when compiling software" OFF)

if(ENABLE_BUILD_CLANG)
    ExternalProject_Add(LLVM
        PREFIX "${CMAKE_BINARY_DIR}"
        #GIT_REPOSITORY "http://llvm.org/git/llvm"
        #GIT_TAG release_37
        DOWNLOAD_COMMAND wget http://llvm.org/releases/3.7.1/llvm-3.7.1.src.tar.xz && tar xvf llvm-3.7.1.src.tar.xz 
        && cp -r "${CMAKE_BINARY_DIR}/src/llvm-3.7.1.src/" ${CMAKE_BINARY_DIR}/src/LLVM
        && wget http://llvm.org/releases/3.7.1/compiler-rt-3.7.1.src.tar.xz && tar xvf compiler-rt-3.7.1.src.tar.xz
        && mv ${CMAKE_BINARY_DIR}/src/compiler-rt-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM/projects/compiler-rt
        && wget http://llvm.org/releases/3.7.1/cfe-3.7.1.src.tar.xz && tar xvf cfe-3.7.1.src.tar.xz
        && mv ${CMAKE_BINARY_DIR}/src/cfe-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM/tools/clang
        CMAKE_ARGS -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${LLVM_INSTALL_PREFIX}
    )

endif()

ExternalProject_Add(OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL "http://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-${OPENMPI_VERSION}.tar.gz"
    CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/src/OpenMPI/configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=${OPENMPI_INSTALL_PREFIX}
)

ExternalProject_Add(Gmsh
    DEPENDS OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL "http://gmsh.info/src/gmsh-${GMSH_VERSION}-source.tgz"
    CMAKE_ARGS -DENABLE_FLTK=${ENABLE_GUI} -DENABLE_MPI=OFF -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=ON -DCMAKE_BUILD_TYPE=release -DENABLE_PETSC=OFF -DCMAKE_INSTALL_PREFIX=${GMSH_INSTALL_PREFIX}
)

ExternalProject_Add(Boost
    DEPENDS OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL "http://sourceforge.net/projects/boost/files/boost/${BOOST_VERSION}/boost_${BOOST_VERSION_US}.tar.gz/download"
    PATCH_COMMAND touch ${CMAKE_BINARY_DIR}/src/Boost/user-config.jam && echo "using mpi ;" >> ${CMAKE_BINARY_DIR}/src/Boost/user-config.jam && echo "" >> ${CMAKE_BINARY_DIR}/src/Boost/user-config.jam
    CONFIGURE_COMMAND cd ${CMAKE_BINARY_DIR}/src/Boost && ./bootstrap.sh
    BUILD_COMMAND cd ${CMAKE_BINARY_DIR}/src/Boost && ./bjam -j16 install --layout=tagged --prefix=${BOOST_INSTALL_PREFIX} --user-config=user-config.jam variant=release threading=single,multi link=static,shared
    INSTALL_COMMAND true
)

ExternalProject_Add(PETSc
    DEPENDS OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz"
    CONFIGURE_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc && PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc ./configure --with-shared-libraries=1  --with-debugging=0  --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3"  --prefix=${PETSC_INSTALL_PREFIX} --download-suitesparse=1  --download-ml  --download-metis  --download-parmetis  --download-blacs  --download-scalapack  --download-fblaslapack  --download-mumps  --download-hypre  --download-ptscotch --download-elemental --download-elemental-shared=1 --with-cxx-dialect=C++11
    BUILD_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc && make PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc PETSC_ARCH=arch-linux2-c-opt all
    INSTALL_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc &&  make PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc PETSC_ARCH=arch-linux2-c-opt install
)

ExternalProject_Add(SLEPc
    DEPENDS PETSc
    PREFIX "${CMAKE_BINARY_DIR}"
    URL "http://slepc.upv.es/download/download.php?filename=slepc-${SLEPC_VERSION}.tar.gz"
    CONFIGURE_COMMAND cd ${CMAKE_BINARY_DIR}/src/SLEPc && SLEPC_DIR=${CMAKE_BINARY_DIR}/src/SLEPc PETSC_DIR=${PETSC_INSTALL_PREFIX} ./configure --prefix=${SLEPC_INSTALL_PREFIX}
    BUILD_COMMAND cd ${CMAKE_BINARY_DIR}/src/SLEPc && make SLEPC_DIR=${CMAKE_BINARY_DIR}/src/SLEPc PETSC_DIR=${PETSC_INSTALL_PREFIX}
    INSTALL_COMMAND cd ${CMAKE_BINARY_DIR}/src/SLEPc && make SLEPC_DIR=${CMAKE_BINARY_DIR}/src/SLEPc PETSC_DIR=${PETSC_INSTALL_PREFIX} install
)

ExternalProject_Add(ParaView
    DEPENDS OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL "http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.4&type=source&os=all&downloadFile=ParaView-v4.4.0-source.tar.gz"
    CMAKE_ARGS -DPARAVIEW_BUILD_QT_GUI=${ENABLE_GUI} -DPARAVIEW_ENABLE_CATALYST=ON -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_USE_MPI=ON -DCMAKE_INSTALL_PREFIX=${PARAVIEW_INSTALL_PREFIX}
)
