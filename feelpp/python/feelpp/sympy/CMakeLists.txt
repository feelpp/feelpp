
find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

feelpp_add_pymodule(sympy SRCS sympy.cpp DESTINATION feelpp/sympy )
file(GLOB_RECURSE SYMPY_PY_FILES
     "${CMAKE_CURRENT_SOURCE_DIR}/*.py"
)

install(FILES ${SYMPY_PY_FILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/sympy)