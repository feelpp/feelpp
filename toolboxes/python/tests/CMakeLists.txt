include(GNUInstallDirs)
include("${CMAKE_CURRENT_LIST_DIR}/../../../feelpp/cmake/modules/FeelppPythonSupport.cmake")

set(_feelpp_pythonpath_entries "${FEELPP_PYTHON_BUILD_DIR}")
list(APPEND _feelpp_pythonpath_entries
  "${CMAKE_BINARY_DIR}/feelpp/python"
  "${CMAKE_BINARY_DIR}/toolboxes/python"
)
list(REMOVE_DUPLICATES _feelpp_pythonpath_entries)

set(_feelpp_ld_library_path_entries)
foreach(_feelpp_ld_dir
    "${CMAKE_BINARY_DIR}/feelpp/feel"
    "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelcore"
    "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelmesh")
  if(EXISTS "${_feelpp_ld_dir}")
    list(APPEND _feelpp_ld_library_path_entries "${_feelpp_ld_dir}")
  endif()
endforeach()
if(FEELPP_TOOLBOXES_ENABLE_CFPDE AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/coefficientformpdes")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/coefficientformpdes")
endif()
if(FEELPP_TOOLBOXES_ENABLE_ELECTRIC AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/electric")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/electric")
endif()
if(FEELPP_TOOLBOXES_ENABLE_HEAT AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heat")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heat")
endif()
if(FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid")
endif()
if(FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermoelectric")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermoelectric")
endif()
if(FEELPP_TOOLBOXES_ENABLE_HEAT AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heatfluid")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heatfluid")
endif()
list(REMOVE_DUPLICATES _feelpp_ld_library_path_entries)

set(PYFILES
  _case_paths.py
  _heat_test_utils.py
  conftest.py
  test_cfpdes.py
  test_collision.py
  test_electric.py
  test_fluid.py
  test_hdg.py
  test_heat.py
  test_heat_thermo2d.py
  test_heat_time_stepping.py
  test_heatfluid.py
  test_interpolation.py
  test_multitoolboxes.py
  test_thermoelectric.py
  test_toolboxes.py
)
feelpp_install_python_tests(COMPONENT toolboxes FILES ${PYFILES})

message(STATUS "[pyfeelpptoolboxes] toolboxes_core: ${TOOLBOXES_CORE}")
message(STATUS "[pyfeelpptoolboxes]      toolboxes: ${TOOLBOXES}")
foreach(item in ${TOOLBOXES} ${TOOLBOXES_CORE})
  set(LDPATH "${FEELPP_PYTHON_BUILD_DIR}/feel/feelmodels/${item}:${LDPATH}")
endforeach()

# Build list of enabled toolboxes for testing
set(_enabled_toolboxes "")
if(FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS)
  list(APPEND _enabled_toolboxes fluid)
endif()
if(FEELPP_TOOLBOXES_ENABLE_CFPDE)
  list(APPEND _enabled_toolboxes cfpdes)
endif()
if(FEELPP_TOOLBOXES_ENABLE_HEAT)
  list(APPEND _enabled_toolboxes heat)
endif()
if(FEELPP_TOOLBOXES_ENABLE_ELECTRIC)
  list(APPEND _enabled_toolboxes electric)
endif()
if(FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC)
  list(APPEND _enabled_toolboxes thermoelectric)
endif()
if(FEELPP_TOOLBOXES_ENABLE_HEAT AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS)
  list(APPEND _enabled_toolboxes heatfluid)
endif()
if(FEELPP_TOOLBOXES_ENABLE_HEAT)
  list(APPEND _enabled_toolboxes interpolation)
endif()

message(STATUS "[pyfeelpptoolboxes] enabled toolboxes for testing: ${_enabled_toolboxes}")

foreach(toolbox IN LISTS _enabled_toolboxes)
  if(toolbox STREQUAL "heat")
    set(_toolbox_test_specs
      "heat|test_heat.py"
      "heat-thermo2d|test_heat_thermo2d.py"
      "heat-time-stepping|test_heat_time_stepping.py"
    )
  else()
    set(_toolbox_test_specs "${toolbox}|test_${toolbox}.py")
  endif()

  # Run pytest on the specific test file to avoid collecting/importing
  # test files for disabled toolboxes
  foreach(_toolbox_test_spec IN LISTS _toolbox_test_specs)
    string(REPLACE "|" ";" _toolbox_test_spec_parts "${_toolbox_test_spec}")
    list(GET _toolbox_test_spec_parts 0 _toolbox_test_name_suffix)
    list(GET _toolbox_test_spec_parts 1 _toolbox_test_file)

    feelpp_add_python_pytest_test(
      NAME feelpp_toolbox_python-tests-${_toolbox_test_name_suffix}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
      LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
      TEST_ARGS -s ${CMAKE_CURRENT_SOURCE_DIR}/${_toolbox_test_file}
    )
    feelpp_add_python_pytest_test(
      NAME feelpp_toolbox_python-tests-mpi-${_toolbox_test_name_suffix}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      MPI_NP 2
      PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
      LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
      TEST_ARGS -s ${CMAKE_CURRENT_SOURCE_DIR}/${_toolbox_test_file}
    )
  endforeach()
endforeach()
