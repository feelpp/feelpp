###  CMakeLists.txt; coding: utf-8 --- 


set(PYFILES conftest.py test_mor.py)
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

foreach(test mor )
  add_test (NAME feelpp_mor_python-tests-${test}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
          ${PYTHON_EXECUTABLE} -m pytest test_${test}.py -s -vv
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_test (NAME feelpp_mor_python-tests-mpi-${test}
    COMMAND ${CMAKE_COMMAND} -E env
          LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
          mpirun -np 2 ${PYTHON_EXECUTABLE} -m pytest test_${test}.py -s -vv
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endforeach()
