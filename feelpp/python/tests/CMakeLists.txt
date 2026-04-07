include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/modules/FeelppPythonSupport.cmake")

set(_feelpp_core_pythonpath_entries "${FEELPP_PYTHON_BUILD_DIR}")
set(_feelpp_core_ld_library_path_entries "${CMAKE_BINARY_DIR}/feelpp/feel")

set(PYFILES
  conftest.py
  test_alg.py
  test_core.py
  test_discr.py
  test_exporter.py
  test_forms.py
  test_integrate.py
  test_measure.py
  test_mesh.py
  test_meshmover.py
  test_operators.py
  test_remesh.py
  test_remotedata.py
  test_timing.py
)
feelpp_install_python_tests(COMPONENT core FILES ${PYFILES})

foreach( test  alg core discr exporter integrate measure mesh meshmover operators remesh remotedata timing )
feelpp_add_python_pytest_test(
  NAME feelpp_qs_python-${test}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  PYTHONPATH_ENTRIES ${_feelpp_core_pythonpath_entries}
  LD_LIBRARY_PATH_ENTRIES ${_feelpp_core_ld_library_path_entries}
  TEST_ARGS -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
)

feelpp_add_python_pytest_test(
  NAME feelpp_qs_python-mpi-${test}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  MPI_NP 2
  PYTHONPATH_ENTRIES ${_feelpp_core_pythonpath_entries}
  LD_LIBRARY_PATH_ENTRIES ${_feelpp_core_ld_library_path_entries}
  TEST_ARGS -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
)
endforeach()
