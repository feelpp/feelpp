include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../ginac)
add_definitions(-DIN_GINAC)

add_executable(viewgar viewgar.cpp)
target_link_libraries(viewgar ginac ${LIBDL_LIBRARIES})
install(TARGETS viewgar RUNTIME DESTINATION "${BIN_INSTALL_DIR}")

if (CMAKE_COMPILER_IS_GNUCC)
	set (CC gcc)
	configure_file (
          "${CMAKE_CURRENT_SOURCE_DIR}/ginac-excompiler.in"
          "${CMAKE_CURRENT_BINARY_DIR}/ginac-excompiler"
	)

	install(PROGRAMS ginac-excompiler DESTINATION "${LIBEXECDIR}")
endif (CMAKE_COMPILER_IS_GNUCC)
