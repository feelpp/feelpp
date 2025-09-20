

set(PYFILES conftest.py test_alg.py test_core.py test_remotedata.py test_discr.py test_exporter.py test_mesh.py test_remesh.py test_meshmover.py)
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

foreach( test  alg core discr exporter integrate measure mesh meshmover operators remesh remotedata timing )
<<<<<<< HEAD
  add_test (NAME feelpp_qs_python-${test}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=:${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}:$ENV{PYTHONPATH} #${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}/..
          ${Python3_EXECUTABLE} -m pytest -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
  )
  #label it so we can do:  ctest -L qs‑python
  set_tests_properties(feelpp_qs_python-${test} PROPERTIES 
    LABELS "qs‑python"
    DEPENDS pyfeelpp

  )

  add_test (NAME feelpp_qs_python-mpi-${test}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=:${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}:$ENV{PYTHONPATH} #${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}/..
          mpirun -np 2 ${Python3_EXECUTABLE} -m pytest -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
  )
  set_tests_properties(feelpp_qs_python-mpi-${test}  PROPERTIES 
    LABELS "qs‑python;mpi"
    DEPENDS pyfeelpp
  )
=======
add_test (NAME feelpp_qs_python-${test}
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=:${CMAKE_BINARY_DIR}/feelpp/feel:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_BINARY_DIR}/python/pyfeelpp:$ENV{PYTHONPATH}
        ${Python3_EXECUTABLE} -m pytest -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test (NAME feelpp_qs_python-mpi-${test}
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=:${CMAKE_BINARY_DIR}/feelpp/feel:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_BINARY_DIR}/python/pyfeelpp:$ENV{PYTHONPATH}
        mpirun -np 2 ${Python3_EXECUTABLE} -m pytest -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
)
>>>>>>> origin/develop
endforeach()