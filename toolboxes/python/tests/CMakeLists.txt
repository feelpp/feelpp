include(GNUInstallDirs)
include("${CMAKE_CURRENT_LIST_DIR}/../../../feelpp/cmake/modules/FeelppPythonSupport.cmake")

function(feelpp_add_python_script_test)
  cmake_parse_arguments(FEELPP_PYSCRIPT
    ""
    "NAME;WORKING_DIRECTORY;MPI_NP;SCRIPT"
    "TEST_ARGS;PYTHONPATH_ENTRIES;LD_LIBRARY_PATH_ENTRIES"
    ${ARGN}
  )
  if(NOT FEELPP_PYSCRIPT_NAME)
    message(FATAL_ERROR "[pyfeelpp] feelpp_add_python_script_test requires NAME")
  endif()
  if(NOT FEELPP_PYSCRIPT_SCRIPT)
    message(FATAL_ERROR "[pyfeelpp] feelpp_add_python_script_test requires SCRIPT")
  endif()

  feelpp_python_make_test_home(_feelpp_python_test_home NAME "${FEELPP_PYSCRIPT_NAME}")
  feelpp_python_join_paths(_feelpp_pythonpath ${FEELPP_PYSCRIPT_PYTHONPATH_ENTRIES})
  feelpp_python_join_paths(_feelpp_python_ld_library_path ${FEELPP_PYSCRIPT_LD_LIBRARY_PATH_ENTRIES})

  set(_feelpp_pyscript_command
    ${CMAKE_COMMAND} -E env
    HOME=${_feelpp_python_test_home}
    PYTHONNOUSERSITE=1
  )

  if(_feelpp_python_ld_library_path)
    list(APPEND _feelpp_pyscript_command
      LD_LIBRARY_PATH=${_feelpp_python_ld_library_path}
    )
  endif()

  if(_feelpp_pythonpath)
    list(APPEND _feelpp_pyscript_command
      PYTHONPATH=${_feelpp_pythonpath}
    )
  endif()

  if(FEELPP_PYSCRIPT_MPI_NP)
    list(APPEND _feelpp_pyscript_command mpirun -np ${FEELPP_PYSCRIPT_MPI_NP})
  endif()

  list(APPEND _feelpp_pyscript_command
    ${Python3_EXECUTABLE}
    ${FEELPP_PYSCRIPT_SCRIPT}
    ${FEELPP_PYSCRIPT_TEST_ARGS}
  )

  if(FEELPP_PYSCRIPT_WORKING_DIRECTORY)
    add_test(
      NAME ${FEELPP_PYSCRIPT_NAME}
      COMMAND ${_feelpp_pyscript_command}
      WORKING_DIRECTORY ${FEELPP_PYSCRIPT_WORKING_DIRECTORY}
    )
  else()
    add_test(NAME ${FEELPP_PYSCRIPT_NAME} COMMAND ${_feelpp_pyscript_command})
  endif()

  if(FEELPP_PYSCRIPT_MPI_NP)
    set_property(TEST ${FEELPP_PYSCRIPT_NAME} PROPERTY PROCESSORS ${FEELPP_PYSCRIPT_MPI_NP})
  endif()
endfunction()

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
if(FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/solid")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/solid")
endif()
if(FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid")
endif()
if(FEELPP_TOOLBOXES_ENABLE_HDG AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/hdg")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/hdg")
endif()
if(FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermoelectric")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/thermoelectric")
endif()
if(FEELPP_TOOLBOXES_ENABLE_HEAT AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heatfluid")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heatfluid")
endif()
list(REMOVE_DUPLICATES _feelpp_ld_library_path_entries)

set(PYFILES
  _cfpdes_test_utils.py
  _cfpdes_case_runner.py
  _case_paths.py
  _heat_test_utils.py
  _solid_test_utils.py
  conftest.py
  test_cfpdes.py
  test_cfpdes_heat_thermal_bridges.py
  test_cfpdes_heat_thermo2d.py
  test_cfpdes_p_laplacian.py
  test_cfpdes_remesh.py
  test_cfpdes_square.py
  test_cfpdes_thermoelectric.py
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
  test_solid.py
  test_solid_nafems_le1.py
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
if(FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS)
  list(APPEND _enabled_toolboxes solid)
endif()
if(FEELPP_TOOLBOXES_ENABLE_HDG)
  list(APPEND _enabled_toolboxes hdg)
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
  if(toolbox STREQUAL "cfpdes")
    set(_cfpdes_test_specs
      "cfpdes|test_cfpdes.py"
      "cfpdes-p-laplacian|test_cfpdes_p_laplacian.py"
      "cfpdes-square|test_cfpdes_square.py"
      "cfpdes-thermoelectric|test_cfpdes_thermoelectric.py"
      "cfpdes-heat-thermal-bridges|test_cfpdes_heat_thermal_bridges.py"
      "cfpdes-heat-thermo2d|test_cfpdes_heat_thermo2d.py"
      "cfpdes-remesh|test_cfpdes_remesh.py"
    )

    foreach(_cfpdes_test_spec IN LISTS _cfpdes_test_specs)
      string(REPLACE "|" ";" _cfpdes_test_spec_parts "${_cfpdes_test_spec}")
      list(GET _cfpdes_test_spec_parts 0 _cfpdes_test_name_suffix)
      list(GET _cfpdes_test_spec_parts 1 _cfpdes_test_file)
      set(_cfpdes_test_args --pytest-file ${CMAKE_CURRENT_SOURCE_DIR}/${_cfpdes_test_file})

      feelpp_add_python_script_test(
        NAME feelpp_toolbox_python-tests-${_cfpdes_test_name_suffix}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/_cfpdes_case_runner.py
        PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
        LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
        TEST_ARGS ${_cfpdes_test_args}
      )
      feelpp_add_python_script_test(
        NAME feelpp_toolbox_python-tests-mpi-${_cfpdes_test_name_suffix}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/_cfpdes_case_runner.py
        MPI_NP 2
        PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
        LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
        TEST_ARGS ${_cfpdes_test_args}
      )
    endforeach()
    continue()
  endif()

  if(toolbox STREQUAL "heat")
    set(_toolbox_test_specs
      "heat|test_heat.py"
      "heat-thermo2d|test_heat_thermo2d.py"
      "heat-time-stepping|test_heat_time_stepping.py"
    )
  elseif(toolbox STREQUAL "solid")
    set(_toolbox_test_specs
      "solid|test_solid.py"
      "solid-nafems-le1|test_solid_nafems_le1.py"
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
