include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../ginac)
add_definitions(-DIN_GINAC)

add_executable(viewgar viewgar.cpp)
target_link_libraries(viewgar ginac)
install(TARGETS viewgar RUNTIME DESTINATION "${BIN_INSTALL_DIR}")

configure_file(ginac-excompiler.in ginac-excompiler)
