
include_directories(${Boost_INCLUDE_DIRS}) 
add_executable(functspace myfunctionspace.cpp myfunctionspace.cfg)
add_library(Funct SHARED myfunctionspace.cpp)
target_link_libraries(Funct ${FEELPP_LIBRARIES})
target_link_libraries(functspace Funct)
target_link_libraries(functspace ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

#add_library(Environment SHARED lib.cpp)
#target_link_libraries(Environment ${FEELPP_LIBRARIES})

#:wqadd_library(Mesh SHARED libMesh.cpp)
#target_link_libraries(Mesh  ${FEELPP_LIBRARIES})

add_library(PyFeelpp SHARED libPyFeelpp.cpp)
target_link_libraries(PyFeelpp ${FEELPP_LIBRARIES})

add_library(PyInteg SHARED libPyInteg.cpp)
target_link_libraries(PyInteg ${FEELPP_LIBRARIES})

add_library(PyLapla SHARED libPyLapla.cpp)
target_link_libraries(PyLapla ${FEELPP_LIBRARIES})

add_library(PyMesh SHARED libPyMesh.cpp)
target_link_libraries(PyMesh ${FEELPP_LIBRARIES})
