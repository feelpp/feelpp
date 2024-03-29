
include(GNUInstallDirs)

set(PYFILES conftest.py test_electric.py test_fluid.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

message(STATUS "[pyfeelpptoolboxes] toolboxes_core: ${TOOLBOXES_CORE}")
message(STATUS "[pyfeelpptoolboxes]      toolboxes: ${TOOLBOXES}")
foreach(item in ${TOOLBOXES} ${TOOLBOXES_CORE})
  set(LDPATH "${CMAKE_BINARY_DIR}/toolboxes/pyfeelpp-toolboxes/feel/feelmodels/${item}:${LDPATH}")
endforeach()

foreach(toolbox fluid cfpdes heat electric thermoelectric heatfluid )
  add_test (NAME feelpp_toolbox_python-tests-${toolbox}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
          ${Python3_EXECUTABLE} -m pytest -k ${toolbox} -s ${CMAKE_CURRENT_SOURCE_DIR}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_test (NAME feelpp_toolbox_python-tests-mpi-${toolbox}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
          mpirun -np 2 ${Python3_EXECUTABLE} -m pytest -k ${toolbox} -s ${CMAKE_CURRENT_SOURCE_DIR}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endforeach()