###  CMakeLists.txt; coding: utf-8 --- 

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
    "${CMAKE_BINARY_DIR}/mor/mor/feel/feelmor"
    "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelcore"
    "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/modelmesh")
  if(EXISTS "${_feelpp_ld_dir}")
    list(APPEND _feelpp_ld_library_path_entries "${_feelpp_ld_dir}")
  endif()
endforeach()
if(FEELPP_TOOLBOXES_ENABLE_HEAT AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heat")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/heat")
endif()
if(FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS AND EXISTS "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid")
  list(APPEND _feelpp_ld_library_path_entries "${CMAKE_BINARY_DIR}/toolboxes/feel/feelmodels/fluid")
endif()
list(REMOVE_DUPLICATES _feelpp_ld_library_path_entries)

set(PYFILES conftest.py test_mor.py test_nirb.py)
feelpp_install_python_tests(COMPONENT mor FILES ${PYFILES})

feelpp_add_python_pytest_test(
  NAME feelpp_mor_python-tests-mor
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
  LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
  TEST_ARGS test_mor.py -s -vv
)
feelpp_add_python_pytest_test(
  NAME feelpp_mor_python-tests-mpi-mor
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  MPI_NP 2
  PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
  LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
  TEST_ARGS test_mor.py -s -vv
)

if(FEELPP_TOOLBOXES_ENABLE_HEAT AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS)
  feelpp_add_python_pytest_test(
    NAME feelpp_mor_python-tests-nirb
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
    LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
    TEST_ARGS test_nirb.py -s -vv -k wogreedy
  )
  feelpp_add_python_pytest_test(
    NAME feelpp_mor_python-tests-mpi-nirb
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    MPI_NP 4
    PYTHONPATH_ENTRIES ${_feelpp_pythonpath_entries}
    LD_LIBRARY_PATH_ENTRIES ${_feelpp_ld_library_path_entries}
    TEST_ARGS test_nirb.py -s -vv -k wogreedy
  )
endif()

# add_test (NAME feelpp_mor_python-tests-greedy-nirb
# COMMAND ${CMAKE_COMMAND} -E env
#       LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
#       PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
#       ${Python3_EXECUTABLE} -m pytest test_nirb.py -s -vv -k egreedy
#       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )
# add_test (NAME feelpp_mor_python-tests-mpi-greedy-nirb
#   COMMAND ${CMAKE_COMMAND} -E env
#         LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
#         PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
#         mpirun -np 2 ${Python3_EXECUTABLE} -m pytest test_nirb.py -s -vv -k egreedy
#         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )
