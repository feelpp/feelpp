

set(PYFILES conftest.py test_alg.py test_core.py test_remotedata.py test_discr.py test_exporter.py test_mesh.py test_remesh.py test_meshmover.py)
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

foreach( test  alg core discr exporter integrate measure mesh meshmover operators remesh remotedata timing )
add_test (NAME feelpp_qs_python-${test}
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=:${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}:$ENV{PYTHONPATH} #${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}/..
        ${Python3_EXECUTABLE} -m pytest -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test (NAME feelpp_qs_python-mpi-${test}
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=:${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}:$ENV{PYTHONPATH} #${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}/..
        mpirun -np 2 ${Python3_EXECUTABLE} -m pytest -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
)
endforeach()