add_custom_target(pyfeelpp)

#foreach(MODULE IN ITEMS feelpp alg mesh discr exporter ts vf models )

macro(feelpp_add_pymodule)
 PARSE_ARGUMENTS(FEELPP_PYMODULE
    "NAME;SRCS"
    ""
    ${ARGN}
    )
  CAR(FEELPP_PYMODULE_NAME ${FEELPP_PYMODULE_DEFAULT_ARGS})
  message(STATUS "[pyfeelpp] add pymodule ${FEELPP_PYMODULE_NAME}")
  pybind11_add_module(_${FEELPP_PYMODULE_NAME}  ${FEELPP_PYMODULE_SRCS}  )
  target_include_directories(_${FEELPP_PYMODULE_NAME} PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
  target_link_libraries( _${FEELPP_PYMODULE_NAME} PUBLIC Feelpp::feelpp )
  add_dependencies(pyfeelpp _${FEELPP_PYMODULE_NAME})
  install(TARGETS _${FEELPP_PYMODULE_NAME} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)
endmacro(feelpp_add_pymodule)

feelpp_add_pymodule(feelpp SRCS feelpp.cpp)
feelpp_add_pymodule(core SRCS core.cpp remotedata.cpp )
feelpp_add_pymodule(alg SRCS alg.cpp)
feelpp_add_pymodule(mesh SRCS mesh.cpp)
feelpp_add_pymodule(discr SRCS discr.cpp)
feelpp_add_pymodule(exporter SRCS exporter.cpp hbf.cpp )
feelpp_add_pymodule(ts SRCS ts.cpp)
feelpp_add_pymodule(vf SRCS vf.cpp)
feelpp_add_pymodule(models SRCS models.cpp)

set(PYFILES __init__.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)