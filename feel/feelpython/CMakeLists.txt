if (FEELPP_ENABLE_PYTHON_WRAPPING AND FEELPP_HAS_BOOST_PYTHON)

include_directories(${Boost_INCLUDE_DIRS})

add_library(PyFeelpp SHARED libPyFeelpp.cpp)
target_link_libraries(PyFeelpp feelpp ${FEELPP_LIBRARIES} ${FEELPP_LINK_LIBRARIES})

#add_library(PyFeelpp2 SHARED libPyFeelpp2.cpp)
#target_link_libraries(PyFeelpp2 feelpp ${FEELPP_LIBRARIES})

add_library(PyInteg SHARED libPyInteg.cpp)
target_link_libraries(PyInteg feelpp ${FEELPP_LIBRARIES} ${FEELPP_LINK_LIBRARIES} PyFeelpp)

add_library(PyLapla SHARED libPyLapla.cpp)
target_link_libraries(PyLapla feelpp ${FEELPP_LIBRARIES} ${FEELPP_LINK_LIBRARIES} PyFeelpp)

#add_library(libtest SHARED libtest.cpp)
#target_link_libraries(libtest feelpp ${FEELPP_LIBRARIES})

endif()
