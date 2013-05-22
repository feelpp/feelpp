# - Try to find OpenTURNS
# Once done this will define
#
#  OpenTURNS_FOUND - system has OT
#  OpenTURNS_INCLUDE_DIR - the OT include directory
#  OpenTURNS_INCLUDE_DIRS - the OT include directory and dependencies include directories
#  OpenTURNS_LIBRARY - Where to find the OT library
#  OpenTURNS_LIBRARIES - Link these to use OT
#  OpenTURNS_WRAPPER_DIR - Wrappers directory
#  OpenTURNS_WRAPPER_DEFINITIONS - Compiler switches required for using OT wrapper
#  OpenTURNS_MODULE_DIR - OT module directory
#  OpenTURNS_MODULE_DEFINITIONS - Compiler switches required for using OT module
#  OpenTURNS_SWIG_INCLUDE_DIR - the OT include directory to swig interface
#
#  Copyright (c) 2009 Mathieu Lapointe <lapointe@phimeca.com>
#  Copyright (c) 2010 Julien Schueller <schueller@phimeca.com>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

include (CheckIncludeFileCXX)
include (FindPackageHandleStandardArgs)

#
# Xml2
#
FIND_PACKAGE(LibXml2)

find_program(OT_CONFIG openturns-config QUIET) 
execute_process( COMMAND ${OT_CONFIG} --version
  OUTPUT_VARIABLE ot_version
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# test if variables are not already in cache
if (NOT (OpenTURNS_INCLUDE_DIR
          AND OpenTURNS_SWIG_INCLUDE_DIR
          AND OpenTURNS_INCLUDE_DIRS
          AND OpenTURNS_LIBRARY
          AND OpenTURNS_BIND_LIBRARY
          AND OpenTURNS_LIBRARIES
          AND OpenTURNS_WRAPPER_DIR
          AND OpenTURNS_MODULE_DIR))

  # set include dir
  if (NOT OpenTURNS_INCLUDE_DIR)
    find_path (OpenTURNS_INCLUDE_DIR
      NAMES
        OT.hxx
      PATHS
        /usr/include
        /usr/local/include
        /opt/local/include
        /sw/include
      PATH_SUFFIXES
        openturns
      DOC
        "OpenTURNS include directory"
    )
  endif ()

  # set swig include dir
  if (NOT OpenTURNS_SWIG_INCLUDE_DIR)
    set(OpenTURNS_SWIG_INCLUDE_DIR "${OpenTURNS_INCLUDE_DIR}/swig")
  endif ()

  # dependencies includes
  if (NOT OpenTURNS_INCLUDE_DIRS)
    set (OpenTURNS_INCLUDE_DIRS ${OpenTURNS_INCLUDE_DIR})
    list (APPEND OpenTURNS_INCLUDE_DIRS ${LIBXML2_INCLUDE_DIR})
    list (APPEND OpenTURNS_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
  endif ()

  # check for library directory
  if (NOT OpenTURNS_LIBRARY)
    find_library (OpenTURNS_LIBRARY
      NAMES
        OT
      PATHS
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /sw/lib
      PATH_SUFFIXES
        openturns
      DOC
        "OpenTURNS library location"
    )
  endif ()
  
  # check for bind library
  if (NOT OpenTURNS_BIND_LIBRARY)
    find_library (OpenTURNS_BIND_LIBRARY
      NAMES
        OTbind
      PATHS
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /sw/lib
      PATH_SUFFIXES
        openturns
      DOC
        "OpenTURNS bind library location"
    )
  endif ()
  
  find_path(wrappers_generic_LIBRARYPATH
      NAMES 
        generic.so
      PATHS 
        /usr/lib/
        /usr/local/lib/
        /opt/local/lib
        /sw/lib/
        PATH_SUFFIXES
        openturns/wrappers
        openturns-${ot_version}/wrappers
      DOC
        "OpenTURNS wrappers generic library location"
  )

  # find dependent libraries
  if (NOT OpenTURNS_LIBRARIES)
    if (NOT OpenTURNS_BIND_LIBRARY)
      set (OpenTURNS_LIBRARIES ${OpenTURNS_LIBRARY} ${LIBXML2_LIBRARIES} ${PYTHON_LIBRARIES})
    else ()
      set (OpenTURNS_LIBRARIES ${OpenTURNS_LIBRARY} ${OpenTURNS_BIND_LIBRARY} ${wrappers_generic_LIBRARYPATH}/generic.so ${LIBXML2_LIBRARIES} ${PYTHON_LIBRARIES})
    endif ()
    list (APPEND OpenTURNS_LIBRARIES ${LIBXML2_LIBRARIES})
    list (APPEND OpenTURNS_LIBRARIES ${PYTHON_LIBRARIES})
  endif ()

  # retrieve path to lib
  get_filename_component (OpenTURNS_LIBRARY_PATH "${OpenTURNS_LIBRARY}" PATH)

  # retrieve install path
  set (OpenTURNS_INSTALL_PATH "${OpenTURNS_LIBRARY_PATH}/../..")

  # find wrappers dir
  if (NOT OpenTURNS_WRAPPER_DIR)
    find_path (OpenTURNS_WRAPPER_DIR
      NAMES
        wrapper.xml wrapper.dtd
      PATHS
        "${OpenTURNS_INSTALL_PATH}"
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /opt/lib
        /opt
      PATH_SUFFIXES
        share/openturns/wrappers
        openturns/wrappers
        openturns-${ot_version}/wrappers
      DOC
        "OpenTURNS wrappers location"
    )
  endif ()

  # set wrapper definitions
  if (NOT OpenTURNS_WRAPPER_DEFINITIONS)
    set(OpenTURNS_WRAPPER_DEFINITIONS)
    check_include_file_cxx (pthread.h FEELPP_HAS_PTHREAD_H)
    if (FEELPP_HAS_PTHREAD_H)
      list (APPEND OpenTURNS_WRAPPER_DEFINITIONS -DFEELPP_HAS_PTHREAD_H)
    endif ()
  endif ()

  # find module directory
  if (NOT OpenTURNS_MODULE_DIR)
    set (OpenTURNS_MODULE_DIR
      "${OpenTURNS_LIBRARY_PATH}/module"
    )
  endif ()

  # set module definitions
  if (NOT OpenTURNS_MODULE_DEFINITIONS)
    set (OpenTURNS_MODULE_DEFINITIONS)

    # check for STDC_HEADERS
    check_include_files (stdlib.h FEELPP_HAS_STDLIB_H)
    check_include_files (stdarg.h FEELPP_HAS_STDARG_H)
    check_include_files (string.h FEELPP_HAS_STRING_H)
    check_include_files (float.h FEELPP_HAS_FLOAT_H)
    check_function_exists (memchr FEELPP_HAS_MEMCHR)
    check_function_exists (free FEELPP_HAS_FREE)
    check_include_files (ctype.h FEELPP_HAS_CTYPE_H)
    if(FEELPP_HAS_STDLIB_H AND FEELPP_HAS_STDARG_H AND FEELPP_HAS_STRING_H AND FEELPP_HAS_FLOAT_H AND FEELPP_HAS_MEMCHR AND FEELPP_HAS_FREE AND FEELPP_HAS_CTYPE_H)
      list (APPEND OpenTURNS_MODULE_DEFINITIONS -DSTDC_HEADERS_H=1)
    else ()
      list (APPEND OpenTURNS_MODULE_DEFINITIONS -DSTDC_HEADERS_H=0)
    endif ()

    # this macro checks a header and defines the corresponding macro
    macro(check_include_files_define_macro header_file)
      # get macro name from header_file
      string(TOUPPER "${header_file}" macro_name)
      string(REGEX REPLACE "[/.]" "_" macro_name ${macro_name})
      set(macro_name FEELPP_HAS_${macro_name})
      # check for header
      check_include_files(${header_file} ${macro_name})
      # define macro
      if(${macro_name})
        list (APPEND OpenTURNS_MODULE_DEFINITIONS -D${macro_name}=1)
      else()
        list (APPEND OpenTURNS_MODULE_DEFINITIONS -D${macro_name}=0)
      endif()
    endmacro()

    # check for some headers
    check_include_files_define_macro(sys/types.h)
    check_include_files_define_macro(sys/stat.h)
    check_include_files_define_macro(stdlib.h)
    check_include_files_define_macro(string.h)
    check_include_files_define_macro(memory.h)
    check_include_files_define_macro(strings.h)
    check_include_files_define_macro(inttypes.h)
    check_include_files_define_macro(stdint.h)
    check_include_files_define_macro(unistd.h)
    check_include_files_define_macro(dlfcn.h)
    check_include_files_define_macro(stdbool.h)
    check_include_files_define_macro(regex.h)

  endif ()

endif ()

# handle the QUIETLY and REQUIRED arguments and set OpenTURNS_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (OPENTURNS DEFAULT_MSG
  OpenTURNS_LIBRARY
  OpenTURNS_BIND_LIBRARY
  OpenTURNS_INCLUDE_DIR
  OpenTURNS_SWIG_INCLUDE_DIR
  OpenTURNS_INCLUDE_DIRS
  OpenTURNS_LIBRARIES
  OpenTURNS_WRAPPER_DIR
  OpenTURNS_MODULE_DIR
)
mark_as_advanced (
  OpenTURNS_LIBRARY
  OpenTURNS_BIND_LIBRARY
  OpenTURNS_INCLUDE_DIR
  OpenTURNS_SWIG_INCLUDE_DIR
  OpenTURNS_INCLUDE_DIRS
  OpenTURNS_LIBRARIES
  OpenTURNS_WRAPPER_DIR
  OpenTURNS_WRAPPER_DEFINITIONS
  OpenTURNS_MODULE_DIR
  OpenTURNS_MODULE_DEFINITIONS
)
