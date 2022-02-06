
include(GNUInstallDirs)

set(PYFILES conftest.py test_mor.py test_reducedbasis.py test_reducedbasis_time.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

message(STATUS "[pyfeelpp-mor] toolboxes_core: ${TOOLBOXES_CORE}") # ??
message(STATUS "[pyfeelpp-mor]      toolboxes: ${TOOLBOXES}")      # ??
foreach(item in ${TOOLBOXES} ${TOOLBOXES_CORE})                    # ??
  set(LDPATH "${CMAKE_BINARY_DIR}/toolboxes/pyfeelpp-toolboxes/feel/feelmodels/${item}:${LDPATH}") # ??
endforeach()


foreach(dim in 2d 3d)
  add_test(NAME pyfeelpp-mor-${dim}
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
        ${PYTHON_EXECUTABLE} -m pytest -k ${dim} -s ${CMAKE_CURRENT_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

  add_test(NAME pyfeelpp-mor-mpi-${dim}
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
        mpirun -np 2 ${PYTHON_EXECUTABLE} -m pytest -k ${dim} -s ${CMAKE_CURRENT_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endforeach()
