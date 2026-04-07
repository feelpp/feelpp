include_guard(GLOBAL)

include(CMakeParseArguments)
include(GNUInstallDirs)

if(NOT DEFINED FEELPP_PYTHON_INSTALL_COMPONENT)
  set(FEELPP_PYTHON_INSTALL_COMPONENT Python)
endif()

function(feelpp_python_make_test_home out_var)
  cmake_parse_arguments(FEELPP_PYHOME
    ""
    "NAME"
    ""
    ${ARGN}
  )
  set(_feelpp_python_test_home "${CMAKE_BINARY_DIR}/python-test-home")
  if(FEELPP_PYHOME_NAME)
    string(REGEX REPLACE "[^A-Za-z0-9_.-]" "_" _feelpp_python_test_home_suffix "${FEELPP_PYHOME_NAME}")
    set(_feelpp_python_test_home "${_feelpp_python_test_home}/${_feelpp_python_test_home_suffix}")
  endif()
  file(MAKE_DIRECTORY "${_feelpp_python_test_home}")
  set(${out_var} "${_feelpp_python_test_home}" PARENT_SCOPE)
endfunction()

function(feelpp_python_filter_existing_paths out_var)
  set(_feelpp_python_paths)
  foreach(_feelpp_python_path IN LISTS ARGN)
    if(_feelpp_python_path AND EXISTS "${_feelpp_python_path}")
      list(APPEND _feelpp_python_paths "${_feelpp_python_path}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES _feelpp_python_paths)
  set(${out_var} "${_feelpp_python_paths}" PARENT_SCOPE)
endfunction()

function(feelpp_python_join_paths out_var)
  set(_feelpp_python_paths ${ARGN})
  list(FILTER _feelpp_python_paths EXCLUDE REGEX "^$")
  list(REMOVE_DUPLICATES _feelpp_python_paths)
  if(_feelpp_python_paths)
    list(JOIN _feelpp_python_paths ":" _feelpp_python_joined_paths)
  else()
    set(_feelpp_python_joined_paths "")
  endif()
  set(${out_var} "${_feelpp_python_joined_paths}" PARENT_SCOPE)
endfunction()

function(feelpp_install_python_tests)
  cmake_parse_arguments(FEELPP_PYTEST_INSTALL
    ""
    "COMPONENT"
    "FILES"
    ${ARGN}
  )
  if(NOT FEELPP_PYTEST_INSTALL_COMPONENT)
    message(FATAL_ERROR "[pyfeelpp] feelpp_install_python_tests requires COMPONENT")
  endif()
  if(FEELPP_PYTEST_INSTALL_FILES)
    install(
      FILES ${FEELPP_PYTEST_INSTALL_FILES}
      DESTINATION ${CMAKE_INSTALL_DATADIR}/feelpp/tests/python/${FEELPP_PYTEST_INSTALL_COMPONENT}
    )
  endif()
endfunction()

function(feelpp_add_python_pytest_test)
  cmake_parse_arguments(FEELPP_PYTEST
    ""
    "NAME;WORKING_DIRECTORY;MPI_NP"
    "TEST_ARGS;PYTHONPATH_ENTRIES;LD_LIBRARY_PATH_ENTRIES"
    ${ARGN}
  )
  if(NOT FEELPP_PYTEST_NAME)
    message(FATAL_ERROR "[pyfeelpp] feelpp_add_python_pytest_test requires NAME")
  endif()

  feelpp_python_make_test_home(_feelpp_python_test_home NAME "${FEELPP_PYTEST_NAME}")
  feelpp_python_join_paths(_feelpp_pythonpath ${FEELPP_PYTEST_PYTHONPATH_ENTRIES})
  feelpp_python_join_paths(_feelpp_python_ld_library_path ${FEELPP_PYTEST_LD_LIBRARY_PATH_ENTRIES})

  set(_feelpp_pytest_command
    ${CMAKE_COMMAND} -E env
    HOME=${_feelpp_python_test_home}
    PYTHONNOUSERSITE=1
  )

  if(_feelpp_python_ld_library_path)
    list(APPEND _feelpp_pytest_command
      LD_LIBRARY_PATH=${_feelpp_python_ld_library_path}
    )
  endif()

  if(_feelpp_pythonpath)
    list(APPEND _feelpp_pytest_command
      PYTHONPATH=${_feelpp_pythonpath}
    )
  endif()

  if(FEELPP_PYTEST_MPI_NP)
    list(APPEND _feelpp_pytest_command mpirun -np ${FEELPP_PYTEST_MPI_NP})
  endif()

  list(APPEND _feelpp_pytest_command ${Python3_EXECUTABLE} -m pytest ${FEELPP_PYTEST_TEST_ARGS})

  if(FEELPP_PYTEST_WORKING_DIRECTORY)
    add_test(
      NAME ${FEELPP_PYTEST_NAME}
      COMMAND ${_feelpp_pytest_command}
      WORKING_DIRECTORY ${FEELPP_PYTEST_WORKING_DIRECTORY}
    )
  else()
    add_test(NAME ${FEELPP_PYTEST_NAME} COMMAND ${_feelpp_pytest_command})
  endif()

  if(FEELPP_PYTEST_MPI_NP)
    set_property(TEST ${FEELPP_PYTEST_NAME} PROPERTY PROCESSORS ${FEELPP_PYTEST_MPI_NP})
  endif()
endfunction()

macro(feelpp_stage_python_files)
  cmake_parse_arguments(FEELPP_STAGE
    ""
    "TARGET;DESTINATION"
    "FILES"
    ${ARGN}
  )
  if(NOT FEELPP_STAGE_TARGET)
    message(FATAL_ERROR "[pyfeelpp] feelpp_stage_python_files requires TARGET")
  endif()
  if(NOT FEELPP_STAGE_DESTINATION)
    message(FATAL_ERROR "[pyfeelpp] feelpp_stage_python_files requires DESTINATION")
  endif()

  if(DEFINED FEELPP_PYTHON_BUILD_DIR)
    set(_FEELPP_STAGE_BUILD_DIR "${FEELPP_PYTHON_BUILD_DIR}")
  else()
    set(_FEELPP_STAGE_BUILD_DIR "${CMAKE_BINARY_DIR}")
  endif()
  get_filename_component(_FEELPP_STAGE_DEST_DIR "${_FEELPP_STAGE_BUILD_DIR}/${FEELPP_STAGE_DESTINATION}" ABSOLUTE)
  file(MAKE_DIRECTORY "${_FEELPP_STAGE_DEST_DIR}")

  # Stage pure Python package files as part of the default build so clean-tree
  # ctest runs do not depend on explicitly building an auxiliary staging target.
  add_custom_target(${FEELPP_STAGE_TARGET} ALL)
  foreach(_FEELPP_STAGE_FILE IN LISTS FEELPP_STAGE_FILES)
    if(IS_ABSOLUTE "${_FEELPP_STAGE_FILE}")
      set(_FEELPP_STAGE_SRC "${_FEELPP_STAGE_FILE}")
    else()
      set(_FEELPP_STAGE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/${_FEELPP_STAGE_FILE}")
    endif()
    get_filename_component(_FEELPP_STAGE_NAME "${_FEELPP_STAGE_FILE}" NAME)
    add_custom_command(
      TARGET ${FEELPP_STAGE_TARGET} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory "${_FEELPP_STAGE_DEST_DIR}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "${_FEELPP_STAGE_SRC}"
              "${_FEELPP_STAGE_DEST_DIR}/${_FEELPP_STAGE_NAME}"
      VERBATIM
    )
  endforeach()
endmacro()

macro(feelpp_add_pymodule)
  set(feelpp_pymodule_name "${ARGV0}")
  cmake_parse_arguments(FEELPP_PYMODULE
    ""
    "DESTINATION;INIT_FILE"
    "SRCS;LINK_LIBRARIES;INCLUDES"
    ${ARGN}
  )
  if(NOT feelpp_pymodule_name)
    message(FATAL_ERROR "[pyfeelpp] feelpp_add_pymodule requires a module name")
  endif()
  if(NOT FEELPP_PYMODULE_DESTINATION)
    message(FATAL_ERROR "[pyfeelpp] feelpp_add_pymodule(${feelpp_pymodule_name}) requires DESTINATION")
  endif()

  message(STATUS "[pyfeelpp] add pymodule ${feelpp_pymodule_name}")
  if(NOT TARGET pybind11::headers)
    if(Python3_EXECUTABLE)
      set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
      set(Python_EXECUTABLE "${Python3_EXECUTABLE}")
    endif()
    if(Python3_INCLUDE_DIRS)
      list(GET Python3_INCLUDE_DIRS 0 PYTHON_INCLUDE_DIR)
      set(PYTHON_INCLUDE_DIRS "${Python3_INCLUDE_DIRS}")
      set(Python_INCLUDE_DIRS "${Python3_INCLUDE_DIRS}")
    endif()
    if(Python3_LIBRARIES)
      list(GET Python3_LIBRARIES 0 PYTHON_LIBRARY)
      set(Python_LIBRARIES "${Python3_LIBRARIES}")
    endif()
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG QUIET)
    if(NOT TARGET pybind11::headers)
      if(TARGET pybind11::pybind11_headers)
        add_library(pybind11::headers ALIAS pybind11::pybind11_headers)
      elseif(TARGET pybind11::pybind11)
        add_library(pybind11::headers ALIAS pybind11::pybind11)
      elseif(TARGET pybind11::module)
        add_library(pybind11::headers ALIAS pybind11::module)
      elseif(TARGET pybind11)
        add_library(pybind11::headers ALIAS pybind11)
      endif()
    endif()
  endif()

  pybind11_add_module(_${feelpp_pymodule_name} ${FEELPP_PYMODULE_SRCS})
  set(_feelpp_python_include_dirs
    ${Python3_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    ${MPI4PY_INCLUDE_DIR}
    ${PETSC4PY_INCLUDE_DIR}
    ${SLEPC4PY_INCLUDE_DIR}
    ${FEELPP_PYMODULE_INCLUDES}
  )
  list(FILTER _feelpp_python_include_dirs EXCLUDE REGEX "^$")
  if(_feelpp_python_include_dirs)
    target_include_directories(
      _${feelpp_pymodule_name}
      PRIVATE
      ${_feelpp_python_include_dirs}
    )
  endif()
  target_link_libraries(
    _${feelpp_pymodule_name}
    PUBLIC
    Feelpp::feelpp
    ${FEELPP_PYMODULE_LINK_LIBRARIES}
  )

  if(DEFINED FEELPP_PYTHON_BUILD_DIR)
    set(_FEELPP_PYTHON_BUILD_DIR "${FEELPP_PYTHON_BUILD_DIR}")
  else()
    set(_FEELPP_PYTHON_BUILD_DIR "${CMAKE_BINARY_DIR}")
  endif()
  get_filename_component(DEST_DIR "${_FEELPP_PYTHON_BUILD_DIR}/${FEELPP_PYMODULE_DESTINATION}" ABSOLUTE)
  file(MAKE_DIRECTORY "${DEST_DIR}")
  set_target_properties(
    _${feelpp_pymodule_name}
    PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY "${DEST_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${DEST_DIR}"
  )
  install(
    TARGETS _${feelpp_pymodule_name}
    DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/${FEELPP_PYMODULE_DESTINATION}
    COMPONENT ${FEELPP_PYTHON_INSTALL_COMPONENT}
  )

  set(_feelpp_python_init_file "")
  if(FEELPP_PYMODULE_INIT_FILE)
    set(_feelpp_python_init_file "${FEELPP_PYMODULE_INIT_FILE}")
  elseif(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/__init__.py")
    set(_feelpp_python_init_file "${CMAKE_CURRENT_BINARY_DIR}/__init__.py")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
    set(_feelpp_python_init_file "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
  endif()

  if(_feelpp_python_init_file)
    add_custom_command(
      TARGET _${feelpp_pymodule_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${DEST_DIR}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${_feelpp_python_init_file}
              ${DEST_DIR}/__init__.py
      VERBATIM
    )
  endif()

  set(_feelpp_python_module_pyfile "")
  if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${feelpp_pymodule_name}.py")
    set(_feelpp_python_module_pyfile "${CMAKE_CURRENT_BINARY_DIR}/${feelpp_pymodule_name}.py")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${feelpp_pymodule_name}.py")
    set(_feelpp_python_module_pyfile "${CMAKE_CURRENT_SOURCE_DIR}/${feelpp_pymodule_name}.py")
  endif()

  if(_feelpp_python_module_pyfile)
    add_custom_command(
      TARGET _${feelpp_pymodule_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${_feelpp_python_module_pyfile}
              ${DEST_DIR}/${feelpp_pymodule_name}.py
      VERBATIM
    )
  endif()
endmacro()
