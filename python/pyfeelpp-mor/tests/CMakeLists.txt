###  CMakeLists.txt; coding: utf-8 --- 


set(PYFILES conftest.py test_mor.py test_nirb.py)
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

add_test (NAME feelpp_mor_python-tests-mor
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
        ${Python3_EXECUTABLE} -m pytest test_mor.py -s -vv
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_test (NAME feelpp_mor_python-tests-mpi-mor
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
        mpirun -np 2 ${Python3_EXECUTABLE} -m pytest test_mor.py -s -vv
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test (NAME feelpp_mor_python-tests-nirb
COMMAND ${CMAKE_COMMAND} -E env
      LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
      PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
      ${Python3_EXECUTABLE} -m pytest test_nirb.py -s -vv -k wogreedy
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_test (NAME feelpp_mor_python-tests-mpi-nirb
COMMAND ${CMAKE_COMMAND} -E env
      LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
      PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
      mpirun -np 2 ${Python3_EXECUTABLE} -m pytest test_nirb.py -s -vv -k wogreedy
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

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
