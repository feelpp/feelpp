

set(PYFILES conftest.py test_alg.py test_core.py test_remotedata.py test_discr.py test_exporter.py test_mesh.py test_remesh.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

add_test (NAME feelpp_qs_python-tests
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=:${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}:$ENV{PYTHONPATH} #${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}/..
        ${PYTHON_EXECUTABLE} -m pytest -vv ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test (NAME feelpp_qs_python-mpi-tests
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=:${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}:$ENV{PYTHONPATH} #${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}/..
        mpirun -np 2 ${PYTHON_EXECUTABLE} -m pytest -vv ${CMAKE_CURRENT_SOURCE_DIR}
)
