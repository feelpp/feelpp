cmake_minimum_required(VERSION 3.19)
project(FeelppContrib NONE)

include(FetchContent)
message(STATUS "${CMAKE_MODULE_PATH}")

#-------------------------------------------------------------------------------
# Generic fetch-or-find helper
#-------------------------------------------------------------------------------
function(find_or_fetch)
    if(ARGC LESS 5)
      message(FATAL_ERROR
        "find_or_fetch(PKG_OPTION PKG_NAME TARGET_NAME GIT_REPO GIT_TAG [VAR=VAL...])\n"
        "  called with ${ARGC} args, but needs at least 5.")
    endif()

    # first five args
    list(GET ARGV 0 PKG_OPTION)
    list(GET ARGV 1 PKG_NAME)
    list(GET ARGV 2 TARGET_NAME)
    list(GET ARGV 3 GIT_REPO)
    list(GET ARGV 4 GIT_TAG)

    # anything beyond 5th is a VAR=VALUE pair
    if(ARGC GREATER 5)
      list(SUBLIST ARGV 5 -1 CONFIG_VARS)
    else()
      set(CONFIG_VARS)
    endif()

    option(FEELPP_USE_EXTERNAL_${PKG_OPTION}
           "Use external ${PKG_NAME}" OFF)

    if(FEELPP_USE_EXTERNAL_${PKG_OPTION})
      message(STATUS "⟳  Using external ${PKG_NAME}")
      find_package(${PKG_NAME} REQUIRED PATHS $ENV{CMAKE_PREFIX_PATH}/lib/cmake/${PKG_NAME})
      message(STATUS "target name: ${TARGET_NAME}")
    else()
      message(STATUS "⇣  Fetching ${PKG_NAME}")
      if(NOT DEFINED FETCHCONTENT_SOURCE_DIR_${PKG_NAME})
        set(FETCHCONTENT_SOURCE_DIR_${PKG_NAME}
            "${CMAKE_CURRENT_SOURCE_DIR}/${PKG_NAME}")
      endif()
      foreach(cfg IN LISTS CONFIG_VARS)
        string(FIND "${cfg}" "=" _eq)
        if(_eq EQUAL -1)
          message(FATAL_ERROR "Bad CONFIG_VARS entry: '${cfg}' (missing '=')")
        endif()
        string(SUBSTRING "${cfg}" 0 ${_eq} _var)
        math(EXPR _ofs "${_eq} + 1")
        string(SUBSTRING "${cfg}" ${_ofs} -1 _val)
        message(STATUS "  - setting ${_var}=${_val}")
        set(${_var} "${_val}" CACHE STRING "" FORCE)
      endforeach()
      FetchContent_Declare(
        ${PKG_NAME}
        GIT_REPOSITORY ${GIT_REPO}
        GIT_TAG        ${GIT_TAG}
      )
      FetchContent_MakeAvailable(${PKG_NAME})
    endif()

    
    if(TARGET ${TARGET_NAME})
      message(STATUS "link with ${TARGET_NAME}")
      target_link_libraries(feelpp_contrib INTERFACE ${TARGET_NAME})
      target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_${PKG_OPTION})
    else()
      message(WARNING " -> TARGET ${TARGET_NAME} does not exist")
    endif()
endfunction()

#----------------------------------------
# Umbrella switch: force all to external
#----------------------------------------
option(FEELPP_USE_EXTERNAL_CONTRIB
       "Force all contrib packages to use external mode"
       OFF)
if(FEELPP_USE_EXTERNAL_CONTRIB)
    message(STATUS "=== Forcing all contrib packages to external ===")
    set(_ALL_PKGS
        GFLAGS GLOG CLN GINAC NANOFLANN FMT NAPP
        METIS PARMETIS NLOPT IPOPT PYBIND11 MONGOCXX
        KWSYS FMI4CPP FMILIB
        TABULATE INDICATORS SIMPLE_WEB_SERVER RANGE_V3
        MATPLOT MMG CPR SPECX EIGENRAND BVH EIGEN3
    )
    foreach(_p IN LISTS _ALL_PKGS)
        set(FEELPP_USE_EXTERNAL_${_p} ON CACHE BOOL "" FORCE)
    endforeach()
endif()

#----------------------------------------
# Interface library for contribs
#----------------------------------------
add_library(feelpp_contrib INTERFACE)
add_library(Feelpp::feelpp_contrib ALIAS feelpp_contrib)

#----------------------------------------
# Fetch or find each third-party dependency
#----------------------------------------
find_or_fetch(NAPP        napp        napp::napp           https://github.com/vincentchabannes/napp.git       v0.3.0)
find_or_fetch(NANOFLANN   nanoflann   nanoflann::nanoflann https://github.com/jlblancoc/nanoflann.git    v1.6.2)
find_or_fetch(FMT         fmt         fmt::fmt             https://github.com/fmtlib/fmt.git               11.0.2)
find_or_fetch(GFLAGS      gflags      gflags::gflags       https://github.com/gflags/gflags.git            v2.2.2)
find_or_fetch(GLOG        glog        glog::glog           https://github.com/google/glog.git              v0.7.1)
find_or_fetch(CLN         cln         cln::cln             https://github.com/feelpp/cln.git               cln_1-3-7)
find_or_fetch(GINAC       ginac       ginac::ginac         https://github.com/feelpp/ginac.git             release-1.8.8)

# Eigen3 via fetch-or-find, with all your original config-vars:
find_or_fetch(
  EIGEN3 Eigen3 Eigen3::Eigen https://gitlab.com/libeigen/eigen/ d6689a15d75e0af62a5f0f6112f6b75c50947eea
#  CONFIG_VARS
    EIGEN_BUILD_DOC=OFF
    BUILD_TESTING=OFF
    EIGEN_LEAVE_TEST_IN_ALL_TARGET=0
    EIGEN_TEST_NOQT=1
)

find_or_fetch(
  MMG mmg mmg::mmg https://github.com/feelpp/mmg.git v7.1.0
#  CONFIG_VARS
    USE_VTK=OFF
    MMG_CI=OFF
    MMG2D_CI=OFF
    MMG3D_CI=OFF
    MMGS_CI=OFF
    BUILD_SHARED_LIBS=ON
    ONLY_VERY_SHORT_TESTS=ON
)

find_or_fetch(
  CPR cpr cpr::cpr https://github.com/feelpp/cpr.git e554da6bde33edc35f1c28bd8197ae50bdf4a4d4
#  CONFIG_VARS
    USE_SYSTEM_CURL=ON
    BUILD_CPR_TESTS=OFF
)

if ( 0 )
find_or_fetch(
  SPECX specx specx::specx https://github.com/feelpp/specx.git eaf74013246099470404ae5b8fd8c7c5b6ad1c5e
#  CONFIG_VARS
    SPECX_JUST_LIB=ON
    SPECX_DO_INSTALL=ON
    SPECX_COMPILE_WITH_MPI=OFF
)
endif( 0 )

find_or_fetch(
  EIGENRAND EigenRand EigenRand::EigenRand https://github.com/feelpp/EigenRand.git 0f25f9ef26d48ec09510d2e8ff07cb12e8f3fa81
#  CONFIG_VARS
    EIGENRAND_BUILD_TEST=OFF
    EIGENRAND_BUILD_BENCHMARK=OFF
)

find_or_fetch(BVH     bvh     Feelpp::feelpp_bvh https://github.com/feelpp/bvh.git 327b67b1d8d798e4dab1686b0dd022a1e0cbcefe)

# Non-fetchable, find-only packages
option(FEELPP_ENABLE_MKL "Enable MKL support" ON)
if(FEELPP_ENABLE_MKL)
  find_package(MKL)
  if(MKL_FOUND)
    target_include_directories(feelpp_contrib INTERFACE ${MKL_INCLUDE_DIRS})
    target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_MKL)
  endif()
endif()



# Metis (find-only)
#option(FEELPP_ENABLE_METIS "Enable Metis support" ON)
#if(FEELPP_ENABLE_METIS)
#  find_package(METIS REQUIRED)
##  if(METIS_FOUND)
##    target_link_libraries(feelpp_contrib INTERFACE metis::metis)
##    target_compile_definitions(feelpp_contrib INTERFACE FEELPP_HAS_METIS)
##  endif()
#endif()
find_or_fetch(METIS     Metis  Metis::metis  https://github.com/numpex/metis.git  v5.3.1-rc.6)


# NLOPT
find_or_fetch(
  NLOPT NLopt nlopt::nlopt https://github.com/stevengj/nlopt.git v2.7.1
 # CONFIG_VARS
    NLOPT_CXX=ON
    NLOPT_FORTRAN=OFF
    BUILD_SHARED_LIBS=ON
    NLOPT_PYTHON=ON
    NLOPT_OCTAVE=OFF
    NLOPT_MATLAB=OFF
    NLOPT_GUILE=OFF
    NLOPT_SWIG=ON
    NLOPT_TESTS=OFF
)

# IPOPT
#find_or_fetch(IPOPT ipopt PkgConfig::IPOPT https://github.com/coin-or/Ipopt.git v3.13.4)

# PyBind11
find_or_fetch(
  PYBIND11 pybind11 pybind11::module https://github.com/pybind/pybind11.git v2.6.2
  CONFIG_VARS
    PYBIND11_INSTALL=ON
)

# MongoCXX
#find_or_fetch(MONGOCXX mongocxx mongocxx::mongocxx https://github.com/mongodb/mongo-cxx-driver.git r3.6.5)

# KWSys
#find_or_fetch(KWSYS kwsys feelpp_kwsys https://github.com/feelpp/kwsys.git master)

# FMI4CPP
find_or_fetch(
  FMI4CPP fmi4cpp fmi4cpp::fmi4cpp https://github.com/feelpp/fmi4cpp.git 645b6a35b1c166ace666c5d387d5a88b074ace18
#  CONFIG_VARS
    USE_SYSTEM_CURL=ON
    BUILD_FMI4CPP_TESTS=ON
)

# FMILIB
#find_or_fetch(FMILIB fmilib fmilib https://github.com/NTNU-IHB/FMILibrary.git v2.7.0)

# Tabulate
find_or_fetch(TABULATE tabulate tabulate::tabulate https://github.com/p‑ranav/tabulate.git v1.5)

# Indicators
find_or_fetch(INDICATORS indicators indicators::indicators https://github.com/p‑ranav/indicators.git v1.9.0)

# Simple Web Server
#find_or_fetch(
#  SIMPLE_WEB_SERVER simple-web-server SimpleWebServer::simple-web-server
#  https://github.com/eidheim/Simple-Web-Server.git v0.5.2
##  CONFIG_VARS
#    BUILD_TESTING=OFF
#)

# Range V3
#find_or_fetch(
#  RANGE_V3 range-v3 range-v3::range-v3 https://github.com/ericniebler/range-v3.git 0.11.0
##  CONFIG_VARS
#    RANGE_V3_INSTALL=ON
#)

# MatPlot++
if (FEELPP_ENABLE_MATPLOT)
  find_or_fetch(MATPLOT Matplot++ matplot::matplot https://github.com/alandefreitas/matplotplusplus.git v1.0)
endif()

message(STATUS "Feelpp contrib libraries configured.")