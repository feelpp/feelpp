###  CMakeLists.txt; coding: utf-8 ---

cmake_minimum_required(VERSION 3.21)
project(toolboxes-multibody)

feelpp_add_pymodule(
  multibody
  SRCS multibody.cpp
  DESTINATION feelpp/toolboxes/multibody
  LINK_LIBRARIES Feelpp::feelpp_toolbox_multibody
)
install(FILES __init__.py
        DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/toolboxes/multibody
        COMPONENT ${FEELPP_PYTHON_INSTALL_COMPONENT})
