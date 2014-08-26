if (FEELPP_ENABLE_PYTHON_WRAPPER) 

include_directories(${Boost_INCLUDE_DIRS}) 
add_executable(myfunctspace myfunctionspace.cpp myfunctionspace.cfg)
add_library(PyFunctSpace SHARED myfunctionspace.cpp)
target_link_libraries(PyFunctSpace feelpp ${FEELPP_LIBRARIES})
target_link_libraries(myfunctspace PyFunctSpace)
target_link_libraries(myfunctspace ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})


add_library(PyFeelpp SHARED libPyFeelpp.cpp)
target_link_libraries(PyFeelpp feelpp ${FEELPP_LIBRARIES})

add_library(PyInteg SHARED libPyInteg.cpp)
target_link_libraries(PyInteg feelpp ${FEELPP_LIBRARIES} PyFeelpp)

add_library(PyLapla SHARED libPyLapla.cpp)
target_link_libraries(PyLapla feelpp ${FEELPP_LIBRARIES} PyFeelpp)

add_library(PyMesh SHARED libPyMesh.cpp)
target_link_libraries(PyMesh feelpp ${FEELPP_LIBRARIES} PyFeelpp)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/mymesh.py" "${CMAKE_CURRENT_BINARY_DIR}/mymesh.py" COPYONLY ) 

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/mylaplacian.py" "${CMAKE_CURRENT_BINARY_DIR}/mylaplacian.py" COPYONLY )
 
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/myintegrals.py" "${CMAKE_CURRENT_BINARY_DIR}/myintegrals.py" COPYONLY )

endif()
