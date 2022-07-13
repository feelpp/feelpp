###  CMakeLists.txt; coding: utf-8 --- 


set(PYFILES conftest.py test_mor.py test_generate_basis.py test_reducedbasis.py)
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

add_test (NAME feelpp_mor_python-tests-mor
COMMAND ${CMAKE_COMMAND} -E env
      LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
      PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
      ${PYTHON_EXECUTABLE} -m pytest test_mor.py -s -vv
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test (NAME feelpp_mor_python-tests-mpi-mor
COMMAND ${CMAKE_COMMAND} -E env
      LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
      PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
      mpirun -np 2 ${PYTHON_EXECUTABLE} -m pytest test_mor.py -s -vv
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# test_generate_basis
foreach (test 2d 3d)

      add_test (NAME feelpp_mor_python-tests-generatebasis-${test}
            COMMAND ${CMAKE_COMMAND} -E env
            LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
            PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
            ${PYTHON_EXECUTABLE} -m pytest test_generate_basis.py -s -vv -k thermal-fin-${test}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )

      add_test (NAME feelpp_mor_python-tests-generatebasis-mpi-${test}
            COMMAND ${CMAKE_COMMAND} -E env
            LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
            PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
            mpirun -np 2 ${PYTHON_EXECUTABLE} -m pytest test_generate_basis.py -s -vv -k thermal-fin-${test}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
endforeach()

# test_reducedbasis
foreach (test 2d 3d)

      add_test (NAME feelpp_mor_python-tests-reducedbasis-${test}
            COMMAND ${CMAKE_COMMAND} -E env
            LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
            PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
            ${PYTHON_EXECUTABLE} -m pytest test_reducedbasis.py -s -vv -k ${test}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )

      add_test (NAME feelpp_mor_python-tests-reducedbasis-mpi-${test}
            COMMAND ${CMAKE_COMMAND} -E env
            LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}
            PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
            mpirun -np 2 ${PYTHON_EXECUTABLE} -m pytest test_reducedbasis.py -s -vv -k ${test}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
endforeach()
