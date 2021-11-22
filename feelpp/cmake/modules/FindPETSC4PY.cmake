# - FindPETSC4PY
# Find PETSC4PY includes and library
# This module defines:
# PETSC4PY_INCLUDE_DIR, where to find PETSC4PY.h, etc.
# PETSC4PY_FOUND

# https://compacc.fnal.gov/projects/repositories/entry/synergia2/CMake/FindPETSC4PY.cmake?rev=c147eafb60728606af4fe7b1b161a660df142e9a

if(NOT PETSC4PY_INCLUDE_DIR)
  execute_process(COMMAND
    ${PYTHON_EXECUTABLE} -c "import petsc4py; print(petsc4py.get_include())"
    OUTPUT_VARIABLE PETSC4PY_INCLUDE_DIR
    RESULT_VARIABLE PETSC4PY_COMMAND_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(PETSC4PY_COMMAND_RESULT)
    message("[feelpp]: PETSC4PY not found ${PETSC4PY_COMMAND_RESULT}")
    set(PETSC4PY_FOUND FALSE)
  else(PETSC4PY_COMMAND_RESULT)
    if (PETSC4PY_INCLUDE_DIR MATCHES "Traceback")
      message("[feelpp]: PETSC4PY matches traceback")
      ## Did not successfully include PETSC4PY
      set(PETSC4PY_FOUND FALSE)
    else (PETSC4PY_INCLUDE_DIR MATCHES "Traceback")
      ## successful
      set(PETSC4PY_FOUND TRUE)
      set(PETSC4PY_INCLUDE_DIR ${PETSC4PY_INCLUDE_DIR} CACHE STRING "PETSC4PY include path")
    endif (PETSC4PY_INCLUDE_DIR MATCHES "Traceback")
  endif(PETSC4PY_COMMAND_RESULT)
else(NOT PETSC4PY_INCLUDE_DIR)
  set(PETSC4PY_FOUND TRUE)
endif(NOT PETSC4PY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC4PY DEFAULT_MSG PETSC4PY_INCLUDE_DIR)
