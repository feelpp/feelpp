if (FEELPP_ENABLE_PYTHON_WRAPPING AND FEELPP_HAS_BOOST_PYTHON) 

include_directories(${Boost_INCLUDE_DIRS}) 
add_executable(myfunctspace myfunctionspace.cpp myfunctionspace.cfg)
add_library(PyFunctSpace SHARED myfunctionspace.cpp)
target_link_libraries(PyFunctSpace feelpp ${FEELPP_LIBRARIES})
target_link_libraries(myfunctspace PyFunctSpace)
target_link_libraries(myfunctspace ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/mymesh.py" "${CMAKE_CURRENT_BINARY_DIR}/mymesh.py" COPYONLY ) 

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/mylaplacian.py" "${CMAKE_CURRENT_BINARY_DIR}/mylaplacian.py" COPYONLY )
 
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/myintegrals.py" "${CMAKE_CURRENT_BINARY_DIR}/myintegrals.py" COPYONLY )

endif()
