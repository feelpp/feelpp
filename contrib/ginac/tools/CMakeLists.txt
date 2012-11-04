include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../ginac)
add_definitions(-DIN_GINAC)

add_executable(viewgar viewgar.cpp)
target_link_libraries(viewgar ginac)
install(TARGETS viewgar RUNTIME DESTINATION "${BIN_INSTALL_DIR}")

configure_file(ginac-excompiler.in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ginac-excompiler)
file(COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ginac-excompiler
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
