add_custom_target(pyfeelpp)

pybind11_add_module(_feelpp  feelpp.cpp  )
target_include_directories(_feelpp PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _feelpp PUBLIC Feelpp::feelpp )

pybind11_add_module(_core  core.cpp remotedata.cpp )
target_include_directories(_core PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _core PUBLIC Feelpp::feelpp )

pybind11_add_module(_mesh  mesh.cpp)
target_include_directories(_mesh PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _mesh PRIVATE Feelpp::feelpp )

pybind11_add_module(_discr  discr.cpp)
target_include_directories(_discr PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _discr PRIVATE  Feelpp::feelpp  )

pybind11_add_module(_exporter  exporter.cpp)
target_include_directories(_exporter PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _exporter PRIVATE  Feelpp::feelpp )

pybind11_add_module(_ts  ts.cpp)
target_include_directories(_ts PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _ts PRIVATE  Feelpp::feelpp  )

pybind11_add_module(_vf  vf.cpp)
target_include_directories(_vf PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _vf PRIVATE Feelpp::feelpp )

pybind11_add_module(_models  models.cpp)
target_include_directories(_models PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
target_link_libraries( _models PRIVATE Feelpp::feelpp )



# pybind11_add_module(pyfeelpp  core.cpp remotedata.cpp  mesh.cpp discr.cpp exporter.cpp ts.cpp vf.cpp models.cpp )
# target_include_directories(pyfeelpp PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR})
# target_link_libraries(pyfeelpp PUBLIC Feelpp::feelpp )
# target_link_libraries(pyfeelpp PRIVATE pybind11::module)
# set_target_properties(pyfeelpp PROPERTIES OUTPUT_NAME feelpp)

add_dependencies(pyfeelpp _feelpp _core _mesh _discr _exporter  _ts _vf _models )

install(TARGETS _feelpp _core _mesh _discr _exporter  _ts _vf _models DESTINATION ${PYTHON_SITE_PACKAGES})

set(PYFILES __init__.py )
install(FILES ${PYFILES} DESTINATION ${PYTHON_SITE_PACKAGES})