
# Check whether a Python module is available by name, and if it is,
# define a variable in the internal cache.
function(_find_python_module_internal module_name version )
  # Check for presence of the module.  Even though we don't use all the
  # variable names set here, assigning them suppresses their output in CMake.

  set(_python_executable "${PYTHON_EXECUTABLE}")
  if ( NOT _python_executable AND Python3_EXECUTABLE )
    set(_python_executable "${Python3_EXECUTABLE}")
  endif()

  if ( NOT _python_executable )
    set(PYTHON_MODULE_${module_name}_FOUND FALSE PARENT_SCOPE)
    return()
  endif()

  execute_process(COMMAND "${_python_executable}" -c "import importlib, re, sys; m=importlib.import_module('${module_name}'); parse=lambda v: tuple(int(x) for x in re.findall(r'\\d+', str(v))[:3]); sys.exit(0 if parse(getattr(m, '__version__', '0')) >= parse('${version}') else 1);"
    RESULT_VARIABLE IMPORT_${module_name}_EXITCODE
    OUTPUT_VARIABLE IMPORT_${module_name}_OUTPUT
    ERROR_VARIABLE IMPORT_${module_name}_ERROR
    )
  
  if(${IMPORT_${module_name}_EXITCODE} EQUAL 0)
    set(PYTHON_MODULE_${module_name}_FOUND TRUE PARENT_SCOPE)
  else()
    set(PYTHON_MODULE_${module_name}_FOUND FALSE PARENT_SCOPE)
  endif()
endfunction()

# Function to simplify checking if a Python module is available
function(find_python_module module_name version result)
  if(NOT PYTHON_MODULE_${module_name}_FOUND)
    _find_python_module_internal(${module_name} ${version})
  endif()
  set(${result} ${PYTHON_MODULE_${module_name}_FOUND} PARENT_SCOPE)
endfunction()
