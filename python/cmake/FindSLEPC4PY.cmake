# - FindSLEPC4PY
# Find SLEPC4PY includes and library
# This module defines:
# SLEPC4PY_INCLUDE_DIR, where to find SLEPC4PY.h, etc.
# SLEPC4PY_FOUND

# https://compacc.fnal.gov/projects/repositories/entry/synergia2/CMake/FindSLEPC4PY.cmake?rev=c147eafb60728606af4fe7b1b161a660df142e9a

if(NOT SLEPC4PY_INCLUDE_DIR)
  execute_process(COMMAND
    ${Python3_EXECUTABLE} -c "import slepc4py; print(slepc4py.get_include())"
    OUTPUT_VARIABLE SLEPC4PY_INCLUDE_DIR
    RESULT_VARIABLE SLEPC4PY_COMMAND_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "SLEPC4PY_INCLUDE_DIR: ${SLEPC4PY_INCLUDE_DIR} command result: ${SLEPC4PY_COMMAND_RESULT}")
  if(SLEPC4PY_COMMAND_RESULT)
    message("[feelpp-python]: SLEPC4PY not found ${SLEPC4PY_COMMAND_RESULT}")
    set(SLEPC4PY_FOUND FALSE)
  else(SLEPC4PY_COMMAND_RESULT)
    if (SLEPC4PY_INCLUDE_DIR MATCHES "Traceback")
      message("[feelpp-python]: SLEPC4PY matches traceback")
      ## Did not successfully include SLEPC4PY
      set(SLEPC4PY_FOUND FALSE)
    else (SLEPC4PY_INCLUDE_DIR MATCHES "Traceback")
      ## successful
      set(SLEPC4PY_FOUND TRUE)
      set(SLEPC4PY_INCLUDE_DIR ${SLEPC4PY_INCLUDE_DIR} CACHE STRING "SLEPC4PY include path")
    endif (SLEPC4PY_INCLUDE_DIR MATCHES "Traceback")
  endif(SLEPC4PY_COMMAND_RESULT)
else(NOT SLEPC4PY_INCLUDE_DIR)
  set(SLEPC4PY_FOUND TRUE)
endif(NOT SLEPC4PY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPC4PY DEFAULT_MSG SLEPC4PY_INCLUDE_DIR)
