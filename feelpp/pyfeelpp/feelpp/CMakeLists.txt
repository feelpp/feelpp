add_custom_target(pyfeelpp)

pybind11_add_module(_core  core.cpp remotedata.cpp  )
target_include_directories(_core PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _core PUBLIC Feelpp::feelpp )
add_dependencies(pyfeelpp _core)
install(TARGETS _core DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)

foreach(MODULE IN ITEMS feelpp alg mesh discr exporter ts vf models )
    pybind11_add_module(_${MODULE}  ${MODULE}.cpp  )
    target_include_directories(_${MODULE} PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
    target_link_libraries( _${MODULE} PUBLIC Feelpp::feelpp )
    add_dependencies(pyfeelpp _${MODULE})
    install(TARGETS _${MODULE} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)
endforeach()

set(PYFILES __init__.py )
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)