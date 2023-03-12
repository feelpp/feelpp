# ##  CMakeLists.txt; coding: utf-8 ---

#set(PYFILES conftest.py)
#install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

add_custom_target(notebooks ALL)

add_dependencies(notebooks venv)

add_test(NAME feelpp-notebook-core
    COMMAND #${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH} CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${PYTHON} -m pytest -s --nbmake ${CMAKE_CURRENT_SOURCE_DIR}/feelpp1.ipynb
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_test(NAME feelpp-notebook-toolbox
    COMMAND # ${CMAKE_COMMAND} -E env        LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:$ENV{LD_LIBRARY_PATH}        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}
    CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${PYTHON} -m pytest -s --nbmake ${CMAKE_CURRENT_SOURCE_DIR}/feelpp2.ipynb
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(feelpp-notebook-core PROPERTIES FIXTURES_REQUIRED notebooks)
set_tests_properties(feelpp-notebook-toolbox PROPERTIES FIXTURES_REQUIRED notebooks)

