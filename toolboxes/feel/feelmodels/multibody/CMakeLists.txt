
feelpp_toolboxes_add_library( multibody
  SRCS multibody.cpp body.cpp
  LINK_LIBRARIES feelpp_modelcore feelpp_modelmesh )

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/multibody COMPONENT Devel)
