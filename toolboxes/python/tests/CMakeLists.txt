
include(GNUInstallDirs)

set(PYFILES conftest.py test_electric.py test_fluid.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

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
if(_enabled_toolboxes)
  list(APPEND _enabled_toolboxes interpolation)
endif()

message(STATUS "[pyfeelpptoolboxes] enabled toolboxes for testing: ${_enabled_toolboxes}")

foreach(toolbox IN LISTS _enabled_toolboxes)
  # Run pytest on the specific test file to avoid collecting/importing
  # test files for disabled toolboxes
  add_test (NAME feelpp_toolbox_python-tests-${toolbox}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
          ${Python3_EXECUTABLE} -m pytest -s ${CMAKE_CURRENT_SOURCE_DIR}/test_${toolbox}.py
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_test (NAME feelpp_toolbox_python-tests-mpi-${toolbox}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
          mpirun -np 2 ${Python3_EXECUTABLE} -m pytest -s ${CMAKE_CURRENT_SOURCE_DIR}/test_${toolbox}.py
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endforeach()
