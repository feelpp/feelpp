
include_directories(${Boost_INCLUDE_DIRS}) 
add_executable(functspace myfunctionspace.cpp myfunctionspace.cfg)
add_library(Funct SHARED myfunctionspace.cpp)
target_link_libraries(Funct ${FEELPP_LIBRARIES})
target_link_libraries(functspace Funct)
target_link_libraries(functspace ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

add_subdirectory(Class)
