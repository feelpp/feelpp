# Internal helper: import MODULE_NAME, check __version__ >= VERSION
function(_find_python_module_internal module_name version)
  # Build a single, properly quoted Python -c command:
  # 1) import pkg_resources to parse versions
  # 2) exit with 0 if module exists and version >= requested
  execute_process(
    COMMAND
      ${PYTHON_EXECUTABLE}
      -c
      "import sys, pkg_resources; \
       dist = pkg_resources.get_distribution('${module_name}'); \
       sys.exit(0 if dist.version and \
                    pkg_resources.parse_version(dist.version) >= \
                    pkg_resources.parse_version('${version}') \
                    else 1)"
    RESULT_VARIABLE IMPORT_${module_name}_EXITCODE
    OUTPUT_VARIABLE IMPORT_${module_name}_OUTPUT
    ERROR_VARIABLE IMPORT_${module_name}_ERROR
  )
  message(STATUS "Checking for Python module ${module_name} >= ${version}: "
                 "exit=${IMPORT_${module_name}_EXITCODE}")
  if(IMPORT_${module_name}_EXITCODE EQUAL 0)
    set(PYTHON_MODULE_${module_name}_FOUND TRUE PARENT_SCOPE)
  else()
    set(PYTHON_MODULE_${module_name}_FOUND FALSE PARENT_SCOPE)
  endif()
endfunction()

# A user‐facing wrapper: only probes once, and returns result in RESULT_VAR
function(find_python_module module_name version result_var)

if (NOT DEFINED PYTHON_EXECUTABLE)
  find_program(PYTHON_EXECUTABLE
    NAMES python3 python
    DOC "Path to the Python executable")
endif()
message(STATUS "Using Python interpreter: ${PYTHON_EXECUTABLE}")
  # If we haven’t already checked, do it now
  if(NOT DEFINED PYTHON_MODULE_${module_name}_FOUND)
    _find_python_module_internal(${module_name} ${version})
  endif()
  # Return the boolean (as CMake TRUE/FALSE)
  set(${result_var} ${PYTHON_MODULE_${module_name}_FOUND} PARENT_SCOPE)
endfunction()