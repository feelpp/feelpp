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
cmake_minimum_required (VERSION 2.8.11)

include(ExternalProject)

include(Deps.cmake)

OPTION(ENABLE_BUILD_CLANG "Enable Clang compilation" OFF)
OPTION(ENABLE_GUI "Enable GUI when compiling software" OFF)

# Let's ensure that we have the correct version of gcc when compiling our applications
set(GNU_INSTALL_DIR "" CACHE PATH "Path where gcc is stored")

set(LOCAL_FILE_PATH "" CACHE PATH "Specify a path where to find packages locally (when internet is unavailable)")
# Allows to add options to BJAM
# Example: if Boost install crashes because it cannot find bzlib.h, you can
# disable bzip2 with -sNO_BZIP2=1
set(BJAM_OPTIONS "" CACHE STRING "Specify options to be passed to bjam")

# Helps with python problems
set(PARAVIEW_OPTIONS "" CACHE STRING "Specify options to be passed to ParaView configure process")

if( NOT GNU_INSTALL_DIR )
    message(FATAL_ERROR "Please specify the installation directory of the gcc distribution")
endif()

if(ENABLE_BUILD_CLANG)
    if(LOCAL_FILE_PATH)
        ExternalProject_Add(LLVM
            PREFIX "${CMAKE_BINARY_DIR}"
            #GIT_REPOSITORY "http://llvm.org/git/llvm"
            #GIT_TAG release_37
            DOWNLOAD_COMMAND 
            rm -rf ${CMAKE_BINARY_DIR}/src/LLVM
            && cp ${LOCAL_FILE_PATH}/llvm-3.7.1.src.tar.xz ${CMAKE_BINARY_DIR}/src 
            && cp ${LOCAL_FILE_PATH}/compiler-rt-3.7.1.src.tar.xz ${CMAKE_BINARY_DIR}/src 
            && cp ${LOCAL_FILE_PATH}/cfe-3.7.1.src.tar.xz ${CMAKE_BINARY_DIR}/src 
            && tar xvf llvm-3.7.1.src.tar.xz && tar xvf compiler-rt-3.7.1.src.tar.xz && tar xvf cfe-3.7.1.src.tar.xz
            && mv ${CMAKE_BINARY_DIR}/src/llvm-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM
            && mv ${CMAKE_BINARY_DIR}/src/compiler-rt-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM/projects/compiler-rt
            && mv ${CMAKE_BINARY_DIR}/src/cfe-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM/tools/clang
            CMAKE_ARGS -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DGCC_INSTALL_PREFIX=${GNU_INSTALL_DIR} -DCMAKE_C_COMPILER=${GNU_INSTALL_DIR}/bin/gcc -DCMAKE_CXX_COMPILER=${GNU_INSTALL_DIR}/bin/g++ -DBUILD_SHARED_LIBS=ON -DLLVM_BUILD_TOOLS=OFF -DCLANG_INCLUDE_DOCS=OFF -DCLANG_INCLUDE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${LLVM_INSTALL_PREFIX}
        )
    else()
        ExternalProject_Add(LLVM
            PREFIX "${CMAKE_BINARY_DIR}"
            #GIT_REPOSITORY "http://llvm.org/git/llvm"
            #GIT_TAG release_37
            DOWNLOAD_COMMAND 
            rm -rf ${CMAKE_BINARY_DIR}/src/LLVM
            && wget http://llvm.org/releases/3.7.1/llvm-3.7.1.src.tar.xz
            && wget http://llvm.org/releases/3.7.1/compiler-rt-3.7.1.src.tar.xz 
            && wget http://llvm.org/releases/3.7.1/cfe-3.7.1.src.tar.xz             
            && tar xvf llvm-3.7.1.src.tar.xz && tar xvf compiler-rt-3.7.1.src.tar.xz && tar xvf cfe-3.7.1.src.tar.xz
            && mv ${CMAKE_BINARY_DIR}/src/llvm-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM
            && mv ${CMAKE_BINARY_DIR}/src/cfe-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM/tools/clang
            && mv ${CMAKE_BINARY_DIR}/src/compiler-rt-3.7.1.src ${CMAKE_BINARY_DIR}/src/LLVM/projects/compiler-rt
            CMAKE_ARGS -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DGCC_INSTALL_PREFIX=${GNU_INSTALL_DIR} -DCMAKE_C_COMPILER=${GNU_INSTALL_DIR}/bin/gcc -DCMAKE_CXX_COMPILER=${GNU_INSTALL_DIR}/bin/g++ -DBUILD_SHARED_LIBS=ON -DLLVM_BUILD_TOOLS=OFF -DCLANG_INCLUDE_DOCS=OFF -DCLANG_INCLUDE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${LLVM_INSTALL_PREFIX}
        )
    endif()
endif()

ExternalProject_Add(OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL ${OPENMPI_URL}
    CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/src/OpenMPI/configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=${OPENMPI_INSTALL_PREFIX}
)

ExternalProject_Add(Gmsh
    DEPENDS OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL ${GMSH_URL}
    CMAKE_ARGS -DENABLE_FLTK=${ENABLE_GUI} -DENABLE_MPI=OFF -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=ON -DCMAKE_BUILD_TYPE=release -DENABLE_PETSC=OFF -DCMAKE_INSTALL_PREFIX=${GMSH_INSTALL_PREFIX}
)

# Ensure that we can pass the argument to the command line
separate_arguments(BJAM_OPTIONS_SEP UNIX_COMMAND "${BJAM_OPTIONS}")

ExternalProject_Add(Boost
    DEPENDS OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL ${BOOST_URL}
    CONFIGURE_COMMAND 
        echo "using mpi : ${OPENMPI_INSTALL_PREFIX}/bin/mpic++ $<SEMICOLON>" > ${CMAKE_BINARY_DIR}/src/Boost/user-config.jam
        && echo "" >> ${CMAKE_BINARY_DIR}/src/Boost/user-config.jam
        && cd ${CMAKE_BINARY_DIR}/src/Boost && ./bootstrap.sh
    BUILD_COMMAND
        cd ${CMAKE_BINARY_DIR}/src/Boost
        && ./bjam -j4 install ${BJAM_OPTIONS_SEP} --layout=tagged --prefix=${BOOST_INSTALL_PREFIX} --user-config=user-config.jam variant=release threading=single,multi link=static,shared
    INSTALL_COMMAND true
)

if(LOCAL_FILE_PATH)
    ExternalProject_Add(PETSc
        DEPENDS OpenMPI
        PREFIX "${CMAKE_BINARY_DIR}"
        URL ${PETSC_URL} 
        CONFIGURE_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc 
            && PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc ./configure 
                --with-shared-libraries=1  --with-debugging=0  --COPTFLAGS="-O3" 
                --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" --prefix=${PETSC_INSTALL_PREFIX} 
                --with-cc=${OPENMPI_INSTALL_PREFIX}/bin/mpicc --with-fc=${OPENMPI_INSTALL_PREFIX}/bin/mpif90 --with-cxx=${OPENMPI_INSTALL_PREFIX}/bin/mpic++
                --download-suitesparse=${LOCAL_FILE_PATH}/SuiteSparse-4.4.3.tar.gz 
                --download-ml=${LOCAL_FILE_PATH}/ml-6.2-p3.tar.gz 
                --download-metis=${LOCAL_FILE_PATH}/metis-5.1.0-p1.tar.gz
                --download-parmetis=${LOCAL_FILE_PATH}/parmetis-4.0.3-p2.tar.gz
                --download-scalapack=${LOCAL_FILE_PATH}/scalapack-2.0.2.tgz 
                --download-fblaslapack=${LOCAL_FILE_PATH}/fblaslapack-3.4.2.tar.gz 
                --download-mumps=${LOCAL_FILE_PATH}/MUMPS_5.0.0-p1.tar.gz 
                --download-hypre=${LOCAL_FILE_PATH}/hypre-2.10.0b-p2.tar.gz 
                --download-ptscotch=${LOCAL_FILE_PATH}/scotch_6.0.3.tar.gz 
                --download-elemental=${LOCAL_FILE_PATH}/Elemental-0.85-p1.tar.gz 
                --download-elemental-shared=1 --with-cxx-dialect=C++11
        BUILD_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc && make PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc PETSC_ARCH=arch-linux2-c-opt all
        INSTALL_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc &&  make PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc PETSC_ARCH=arch-linux2-c-opt install
    )
else()
    ExternalProject_Add(PETSc
        DEPENDS OpenMPI
        PREFIX "${CMAKE_BINARY_DIR}"
        URL ${PETSC_URL} 
        CONFIGURE_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc 
            && PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc ./configure 
                --with-shared-libraries=1  --with-debugging=0  --COPTFLAGS="-O3" 
                --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" --prefix=${PETSC_INSTALL_PREFIX} 
                --with-cc=${OPENMPI_INSTALL_PREFIX}/bin/mpicc --with-fc=${OPENMPI_INSTALL_PREFIX}/bin/mpif90 --with-cxx=${OPENMPI_INSTALL_PREFIX}/bin/mpic++
                --download-suitesparse=1 --download-ml --download-metis --download-parmetis 
                --download-blacs --download-scalapack=1 --download-fblaslapack --download-mumps 
                --download-hypre=1 --download-ptscotch --download-elemental 
                --download-elemental-shared=1 --with-cxx-dialect=C++11
        BUILD_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc && make PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc PETSC_ARCH=arch-linux2-c-opt all
        INSTALL_COMMAND cd ${CMAKE_BINARY_DIR}/src/PETSc &&  make PETSC_DIR=${CMAKE_BINARY_DIR}/src/PETSc PETSC_ARCH=arch-linux2-c-opt install
    )
endif()

ExternalProject_Add(SLEPc
    DEPENDS PETSc
    PREFIX "${CMAKE_BINARY_DIR}"
    URL ${SLEPC_URL}
    CONFIGURE_COMMAND cd ${CMAKE_BINARY_DIR}/src/SLEPc && SLEPC_DIR=${CMAKE_BINARY_DIR}/src/SLEPc PETSC_DIR=${PETSC_INSTALL_PREFIX} ./configure --prefix=${SLEPC_INSTALL_PREFIX}
    BUILD_COMMAND cd ${CMAKE_BINARY_DIR}/src/SLEPc && make SLEPC_DIR=${CMAKE_BINARY_DIR}/src/SLEPc PETSC_DIR=${PETSC_INSTALL_PREFIX}
    INSTALL_COMMAND cd ${CMAKE_BINARY_DIR}/src/SLEPc && make SLEPC_DIR=${CMAKE_BINARY_DIR}/src/SLEPc PETSC_DIR=${PETSC_INSTALL_PREFIX} install
)

# Ensure that we can pass the argument to the command line
separate_arguments(PARAVIEW_OPTIONS_SEP UNIX_COMMAND "${PARAVIEW_OPTIONS}")

ExternalProject_Add(ParaView
    DEPENDS OpenMPI
    PREFIX "${CMAKE_BINARY_DIR}"
    URL ${PARAVIEW_URL}
    CMAKE_ARGS -DBUILD_TESTING=OFF -DPARAVIEW_BUILD_QT_GUI=${ENABLE_GUI} -DPARAVIEW_ENABLE_CATALYST=ON -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_USE_MPI=ON -DCMAKE_INSTALL_PREFIX=${PARAVIEW_INSTALL_PREFIX} -DMPIEXEC=${OPENMPI_INSTALL_PREFIX}/bin/mpiexec ${PARAVIEW_OPTIONS_SEP}
)

# Configure file for downloading packages
# Useful when on an other machine internet is not available and we want to use LOCAL_FILE_PATH
configure_file(get_distrib.sh.in get_distrib.sh @ONLY@)

# Create a profile module for the packages
configure_file(feelpp.profile.in feelpp.profile @ONLY@)

# Create a shell script to source for loading the built applications
configure_file(setup_env.sh.in setup_env.sh @ONLY@)
